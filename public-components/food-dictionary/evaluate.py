"""Compute mapping metrics only on actual reviewed decisions; never invent gold."""
import json
from pathlib import Path
import argparse


def evaluate(rows,k=10):
    reviewed=[r for r in rows if r.get('reviewer') and isinstance(r.get('gold_equivalent'),bool)]
    positive=[r for r in reviewed if r['gold_equivalent']]
    predictions=[r for r in reviewed if r.get('prediction_equivalent') is True]
    decided=[r for r in reviewed if isinstance(r.get('prediction_equivalent'),bool)]
    true_positive=sum(r['gold_equivalent'] for r in predictions)
    false_positive=len(predictions)-true_positive
    ratio=lambda a,b:a/b if b else None
    return {'reviewed_rows':len(reviewed),'gold_positive_rows':len(positive),'accepted_predictions':len(predictions),
            'precision':ratio(true_positive,len(predictions)),'recall':ratio(true_positive,len(positive)),
            'false_merge_rate_among_accepted':ratio(false_positive,len(predictions)),
            'retrieval_recall_at_k':ratio(sum(r.get('object_id') in r.get('retrieved_object_ids',[])[:k] for r in positive),len(positive)),
            'decision_coverage':ratio(len(decided),len(reviewed)),'abstention_rate':ratio(len(reviewed)-len(decided),len(reviewed)),
            'k':k,'human_gold_available':bool(reviewed)}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('review_set',type=Path);a=p.parse_args()
    print(json.dumps(evaluate([json.loads(line) for line in a.review_set.read_text().splitlines()]),indent=2))
