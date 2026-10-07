"""Source-offset gene mentions from a locally selected token classifier.

Mentions are entity predictions. They do not identify canonical genes or imply
compound relationships, dietary effects or article-independent accuracy.
"""
from __future__ import annotations
import re
from pathlib import Path
from .token_evaluation import encode_word_chunks,reconstruct_word_predictions,spans

LABELS=['O','B-GENE','I-GENE']


def native_word_offsets(text):
    # Preserve exact Unicode-codepoint offsets; punctuation remains its own word.
    return [{'word':m.group(),'start':m.start(),'end':m.end()}
            for m in re.finditer(r"\w+(?:[-']\w+)*|[^\w\s]",text,flags=re.UNICODE)]


def mentions_from_labels(text,tokens,labels):
    if len(tokens)!=len(labels):raise ValueError('Source-offset/token label mismatch')
    for token in tokens:
        if text[token['start']:token['end']]!=token['word']:raise ValueError('Token offsets do not match the source text')
    mentions=[]
    for start,end,kind in sorted(spans(labels)):
        if kind!='GENE':raise ValueError('The selected model emits gene entities only')
        lo,hi=tokens[start]['start'],tokens[end-1]['end']
        mentions.append({'start':lo,'end':hi,'text':text[lo:hi],'type':'GENE',
                         'offset_unit':'Unicode codepoints','origin':'BC2GM-finetuned BiomedBERT prediction'})
    return mentions


class NeuralGeneMentionProvider:
    def __init__(self,model_directory,device='cpu'):
        import torch
        from transformers import AutoTokenizer,AutoModelForTokenClassification
        path=Path(model_directory)
        self.tokenizer=AutoTokenizer.from_pretrained(path,local_files_only=True)
        self.model=AutoModelForTokenClassification.from_pretrained(path,local_files_only=True).to(device).eval()
        self.device=device;self.torch=torch
        actual={int(i):v for i,v in self.model.config.id2label.items()}
        if actual!=dict(enumerate(LABELS)):raise ValueError('Selected model label mapping is incompatible')

    def extract(self,text):
        tokens=native_word_offsets(text)
        if not tokens:return []
        rows=[{'words':[r['word'] for r in tokens]}]
        pieces=encode_word_chunks(self.tokenizer,rows,LABELS,False,max_tokens=254);predictions=[]
        with self.torch.inference_mode():
            for piece in pieces:
                inputs={k:self.torch.tensor([v],device=self.device) for k,v in piece['features'].items() if k!='labels'}
                predictions.append(self.model(**inputs).logits.argmax(-1).cpu()[0].tolist())
        labels=reconstruct_word_predictions(rows,pieces,predictions,LABELS)[0]
        return mentions_from_labels(text,tokens,labels)
