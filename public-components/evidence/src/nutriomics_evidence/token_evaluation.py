"""Exact native-word spans and checked wordpiece/chunk reconstruction.

An orphan I tag begins a new span, consistent with IOB2 repair. Exact type and
both boundaries are required. This is native-word evaluation, not document or
character-offset independence.
"""
from collections import Counter


def spans(labels):
    result=set();start=None;kind=None
    for index,label in enumerate(list(labels)+['O']):
        if label != 'O' and (not isinstance(label,str) or not label.startswith(('B-','I-')) or len(label)<=2):
            raise ValueError(f'Invalid IOB entity label: {label!r}')
        tag,_,name=label.partition('-')
        if start is not None and (tag!='I' or name!=kind):result.add((start,index,kind));start=None
        if tag=='B' or (tag=='I' and start is None):start=index;kind=name
    return result


def span_metrics(gold,predicted):
    if len(gold)!=len(predicted) or any(len(a)!=len(b) for a,b in zip(gold,predicted)):raise ValueError('Token/source alignment mismatch')
    tp=fp=fn=0
    for a,b in zip(gold,predicted):
        observed,estimated=spans(a),spans(b);tp+=len(observed&estimated);fp+=len(estimated-observed);fn+=len(observed-estimated)
    precision=tp/(tp+fp) if tp+fp else 0.;recall=tp/(tp+fn) if tp+fn else 0.
    return {'precision':precision,'recall':recall,'f1':2*precision*recall/(precision+recall) if precision+recall else 0.,'true_positive_spans':tp,'false_positive_spans':fp,'false_negative_spans':fn}


def encode_word_chunks(tokenizer,rows,labels,use_gold,max_tokens=254):
    """Split at complete native words and retain global indices, with no loss."""
    if not tokenizer.is_fast: raise ValueError('Fast tokenizer word_ids are required')
    special=tokenizer.num_special_tokens_to_add(pair=False)
    if max_tokens<=special: raise ValueError('No encoding budget remains for native words')
    pieces=[]
    for sid,row in enumerate(rows):
        words=row['words']
        if not words or any(not isinstance(w,str) or not w for w in words): raise ValueError('Empty native sentence or word')
        full=tokenizer(words,is_split_into_words=True,add_special_tokens=False,truncation=False)
        counts=Counter(full.word_ids())
        if any(counts[index]==0 for index in range(len(words))): raise ValueError('Tokenizer discarded a native word')
        start=0
        while start<len(words):
            end=start;length=special
            while end<len(words) and length+counts[end]<=max_tokens: length+=counts[end];end+=1
            if end==start: raise ValueError('A native word exceeds the encoding budget')
            batch=tokenizer(words[start:end],is_split_into_words=True,truncation=False)
            if len(batch['input_ids'])>max_tokens: raise ValueError('Measured chunk exceeded its declared encoding budget')
            mapping=[None if index is None else start+index for index in batch.word_ids()]
            if {w for w in mapping if w is not None}!=set(range(start,end)): raise ValueError('Chunk lost native word coverage')
            seen=set();native=[]
            for word in mapping:
                if word is None or word in seen or not use_gold: native.append(-100)
                else: native.append(labels.index(row['labels'][word]))
                if word is not None: seen.add(word)
            pieces.append({'features':dict(batch)|{'labels':native},'sentence':sid,'mapping':mapping,'start':start,'end':end})
            start=end
    return pieces


def reconstruct_word_predictions(rows,pieces,tag_sequences,labels):
    """Use first wordpiece only; reject silent drops, overlap and short outputs."""
    if len(pieces)!=len(tag_sequences): raise ValueError('Chunk prediction alignment mismatch')
    result=[[None]*len(r['words']) for r in rows]
    for piece,tags in zip(pieces,tag_sequences):
        if len(tags)<len(piece['mapping']): raise ValueError('Truncated chunk prediction')
        seen=set()
        for word,tag in zip(piece['mapping'],tags):
            if word is None or word in seen: continue
            if not isinstance(tag,int) or not 0<=tag<len(labels): raise ValueError('Unknown neural entity tag')
            if result[piece['sentence']][word] is not None: raise ValueError('Overlapping prediction for a native word')
            result[piece['sentence']][word]=labels[tag];seen.add(word)
    if any(label is None for sentence in result for label in sentence): raise ValueError('Unpredicted native word')
    return result
