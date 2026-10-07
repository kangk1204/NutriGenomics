"""Constrained image/text matching, followed by explicit human food/portion confirmation."""
import hashlib
import io
import json
import math
from pathlib import Path
import threading
from PIL import Image, ImageOps, UnidentifiedImageError

LABELS = ('apple','banana','orange','strawberry','grapes','avocado','tomato','carrot','broccoli','spinach','potato','sweet potato','rice','bread','oatmeal','pasta','pizza','chicken','beef','salmon','egg','milk','yogurt','cheese','tofu','coffee','tea','chocolate','almonds','walnuts')


def decode_image(raw):
    if not raw or len(raw)>8*1024*1024:raise ValueError('Image must be between 1 byte and 8 MiB')
    try:
        with Image.open(io.BytesIO(raw)) as image:
            if image.width*image.height>20_000_000:raise ValueError('Image exceeds 20 million pixels')
            image.load();return ImageOps.exif_transpose(image).convert('RGB')
    except (UnidentifiedImageError,OSError,Image.DecompressionBombError) as exc:
        raise ValueError('Image is not a decodable PNG/JPEG/WebP') from exc


class FoodImageMatcher:
    def __init__(self,root,model_dir,food_registry):
        self.root=Path(root).resolve();self.model_dir=Path(model_dir).resolve();self.food_registry=food_registry
        if not self.model_dir.is_relative_to(self.root):raise ValueError('Image model must remain inside research root')
        self.lock=threading.Lock();self.model=None

    def predict(self,raw):
        image=decode_image(raw)
        with self.lock:
            if not self.model_dir.is_dir():raise FileNotFoundError('Pinned image candidate model is not installed; use food search')
            import torch
            from transformers import CLIPModel,CLIPProcessor
            torch.set_num_threads(2)
            if self.model is None:
                receipt_path=self.model_dir.parent/'clip_receipt.json'
                if not receipt_path.is_file():raise FileNotFoundError('Pinned image model receipt is not installed')
                receipt=json.loads(receipt_path.read_text())
                for name,expected in receipt['files'].items():
                    path=(self.model_dir/name).resolve()
                    if not path.is_relative_to(self.model_dir) or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:raise ValueError('Image model/tokenizer checksum mismatch')
                model=CLIPModel.from_pretrained(self.model_dir,local_files_only=True).eval()
                processor=CLIPProcessor.from_pretrained(self.model_dir,local_files_only=True)
                self.model,self.processor=model,processor
            inputs=self.processor(text=['a photo of '+v for v in LABELS],images=image,return_tensors='pt',padding=True)
            with torch.inference_mode():score=self.model(**inputs).logits_per_image[0].softmax(0).tolist()
            if len(score)!=len(LABELS) or any(not isinstance(v,(int,float)) or not math.isfinite(v) or not 0<=v<=1 for v in score):
                raise ValueError('Image model returned invalid relative matching scores')
        top=sorted(zip(LABELS,score),key=lambda x:x[1],reverse=True)[:5]
        items=[]
        for term,value in top:
            # This returns source foods for confirmation, never a measured portion or efficacy.
            for food in self.food_registry.search(term,3)['items']:
                items.append(food|{'image_candidate':term,'relative_matching_score':value})
        return {'items':items,'candidate_terms':[{'term':v,'relative_matching_score':s} for v,s in top],
                'model':'CLIP ViT-B/32','input_sha256':hashlib.sha256(raw).hexdigest(),
                'requires_food_confirmation':True,'requires_portion_confirmation':True,
                'score_interpretation':'Relative match within a fixed 30-food vocabulary; accuracy is not calibrated',
                'clinical_effect_inferred':False,'image_estimated_portion':False}
