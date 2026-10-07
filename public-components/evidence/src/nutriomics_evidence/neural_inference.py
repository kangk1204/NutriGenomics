"""Frozen chemical–gene relation-presence inference with verified source mention offsets."""
import json
import hashlib
from pathlib import Path
import threading
from .neural_data import mark_context,sha


class EvidencePredictor:
    def __init__(self,model_dir,receipt):
        self.model_dir=Path(model_dir);self.receipt=Path(receipt);self.model=None;self.lock=threading.Lock()
        self.receipt_sha256=None;self.frozen_receipt=None

    def predict(self,text,chemical_start,chemical_end,gene_start,gene_end,source_id):
        if not isinstance(text,str) or not 1<=len(text)<=50_000:raise ValueError('Source text must be 1–50,000 characters')
        chemical={'start':chemical_start,'end':chemical_end,'text':text[chemical_start:chemical_end]}
        gene={'start':gene_start,'end':gene_end,'text':text[gene_start:gene_end]}
        marked=mark_context(text,chemical,gene)
        with self.lock:
            if not self.receipt.is_file():raise FileNotFoundError('Frozen evidence model receipt is not installed')
            receipt_bytes=self.receipt.read_bytes()
            receipt_hash=hashlib.sha256(receipt_bytes).hexdigest()
            if self.model is not None and receipt_hash != self.receipt_sha256:
                raise ValueError('Frozen evidence receipt changed after model initialization')
            receipt=json.loads(receipt_bytes)
            if self.model is None:
                if sha(self.model_dir/'model.safetensors')!=receipt['model_sha256']:raise ValueError('Evidence model checksum mismatch')
                for name,expected in receipt.get('model_files',{}).items():
                    path=(self.model_dir/name).resolve()
                    if not path.is_relative_to(self.model_dir.resolve()) or sha(path)!=expected:raise ValueError('Evidence encoder/tokenizer checksum mismatch')
                import torch
                from transformers import AutoTokenizer,AutoModelForSequenceClassification
                torch.set_num_threads(2)
                tokenizer=AutoTokenizer.from_pretrained(self.model_dir,local_files_only=True)
                model=AutoModelForSequenceClassification.from_pretrained(self.model_dir,local_files_only=True).eval()
                self.tokenizer,self.model=tokenizer,model
                self.receipt_sha256,self.frozen_receipt=receipt_hash,receipt
            receipt=self.frozen_receipt
            import torch
            inputs=self.tokenizer(marked,truncation=True,max_length=384,return_tensors='pt')
            with torch.inference_mode():
                probabilities=torch.sigmoid(self.model(**inputs).logits.float())[0]
                score=float(1-torch.prod(1-probabilities))
        return {'source_id':source_id,'source_text_sha256':__import__('hashlib').sha256(text.encode()).hexdigest(),
                'chemical_mention':chemical,'gene_mention':gene,'score':score,
                'predicted_relation':score>=receipt['threshold'],'threshold':receipt['threshold'],
                'model_sha256':receipt['model_sha256'],'protocol_id':receipt['protocol_id'],
                'entity_input':'user-confirmed source offsets',
                'endpoint':'canonical chemical–gene literature relation presence',
                'clinical_efficacy_evaluated':False,'expert_review_required':True}
