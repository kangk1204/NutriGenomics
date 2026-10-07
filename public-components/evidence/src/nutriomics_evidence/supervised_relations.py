"""Source-grounded, CPU relation-extraction comparator with supplied mentions.

This is a supervised text model, not an LLM, NER model or food-domain gold.
Every candidate retains literal source offsets and chemical/gene entity IDs.
"""
from __future__ import annotations
from collections import defaultdict
import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import SGDClassifier
from sklearn.multiclass import OneVsRestClassifier
from .extraction import RELATIONS, sentence_spans

LABELS = tuple(sorted(RELATIONS))


def candidates(document, max_span=800):
    text = document["text"]
    spans = list(sentence_spans(text))
    chemicals = [e for e in document["entities"] if e["type"] == "CHEMICAL"]
    genes = [e for e in document["entities"] if e["type"].startswith("GENE")]
    output = []
    for chemical in chemicals:
        for gene in genes:
            left, right = min(chemical["start"], gene["start"]), max(chemical["end"], gene["end"])
            if right - left > max_span:
                continue
            start = next((a for a, b in spans if a <= left < b), left)
            end = next((b for a, b in spans if a < right <= b), right)
            if end - start > max_span:
                start, end = left, right
            replacements = []
            for entity in document["entities"]:
                if start <= entity["start"] and entity["end"] <= end:
                    marker = "SELECTED_CHEMICAL" if entity["id"] == chemical["id"] else "SELECTED_GENE" if entity["id"] == gene["id"] else "OTHER_CHEMICAL" if entity["type"] == "CHEMICAL" else "OTHER_GENE"
                    replacements.append((entity["start"] - start, entity["end"] - start, " " + marker + " "))
            # Mention inventories with overlapping ranges are unsuitable for
            # this simple masked lexical baseline; retain the literal window.
            replacements.sort()
            overlap = any(a[1] > b[0] for a, b in zip(replacements, replacements[1:]))
            masked = text[start:end]
            if not overlap:
                for a, b, marker in reversed(replacements):
                    masked = masked[:a] + marker + masked[b:]
            output.append({"pmid": document["pmid"], "arg1": chemical["id"], "arg2": gene["id"],
                           "evidence_start": start, "evidence_end": end, "evidence_text": text[start:end],
                           "features": masked, "overlapping_mentions": overlap})
    return output


def targets(rows, relations):
    lookup = defaultdict(set)
    for pmid, label, chemical, gene in relations:
        lookup[(pmid, chemical, gene)].add(label)
    y = np.zeros((len(rows), len(LABELS)), dtype=np.int8)
    for i, row in enumerate(rows):
        for label in lookup[(row["pmid"], row["arg1"], row["arg2"])]:
            y[i, LABELS.index(label)] = 1
    return y


def fit(rows, y, seed=20261003):
    vectorizer = TfidfVectorizer(ngram_range=(1, 2), min_df=2, max_features=150000,
                               sublinear_tf=True, dtype=np.float32)
    X = vectorizer.fit_transform(row["features"] for row in rows)
    model = OneVsRestClassifier(SGDClassifier(loss="log_loss", alpha=1e-5, max_iter=1000,
                                             tol=1e-4, random_state=seed, average=True), n_jobs=1)
    model.fit(X, y)
    return vectorizer, model


def predict(rows, vectorizer, model):
    return model.predict_proba(vectorizer.transform(row["features"] for row in rows))


def export(rows, probabilities, threshold):
    result = []
    if probabilities.shape != (len(rows), len(LABELS)):
        raise ValueError("Candidate/label probability alignment failed")
    if not np.isfinite(probabilities).all():
        raise ValueError("Nonfinite relation probabilities")
    for row, values in zip(rows, probabilities):
        for j in np.flatnonzero(values >= threshold):
            result.append({key: row[key] for key in ("pmid", "arg1", "arg2", "evidence_start", "evidence_end", "evidence_text")}
                          | {"relation": LABELS[j], "probability": float(values[j]),
                             "method": "supervised_tfidf_sgd", "tier": "extracted_requires_review",
                             "binding_affinity": None, "semantic_entailment_validated": False})
    return result


def tuple_scores(gold, rows):
    prediction = {(r["pmid"], r["relation"], r["arg1"], r["arg2"]) for r in rows}
    tp, fp, fn = len(gold & prediction), len(prediction - gold), len(gold - prediction)
    precision, recall = tp / (tp + fp) if tp + fp else 0., tp / (tp + fn) if tp + fn else 0.
    return {"tp": tp, "fp": fp, "fn": fn, "precision": precision, "recall": recall,
            "f1": 2 * precision * recall / (precision + recall) if precision + recall else 0.}
