"""V36: real TEXT encoding counts, not API billing or reasoning-token estimates."""
from pathlib import Path
import os,sys,json

ENCODINGS=('o200k_base','cl100k_base')


def encoder(name='o200k_base'):
    if name not in ENCODINGS:raise ValueError('Choose an explicitly measured text encoding')
    try:import tiktoken
    except ImportError:
        local=Path(__file__).resolve().parents[2]/'work'/'token-deps'
        if local.is_dir():sys.path.insert(0,str(local))
        try:import tiktoken
        except ImportError:raise RuntimeError('Text token counts unavailable: install tiktoken==0.11.0; character estimates are not substituted')
    os.environ.setdefault('TIKTOKEN_CACHE_DIR',str(Path(__file__).resolve().parent/'token_cache'))
    return tiktoken.get_encoding(name)


def wire(value):return json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False)
def count(text,name='o200k_base'):return len(encoder(name).encode(text,disallowed_special=()))
def counts(text):return {name:count(text,name) for name in ENCODINGS}


def conversation(events,instructions='',schema=''):
    """All provided request/response/context text counted, excluding unknown API framing."""
    texts=[instructions,schema]+[event['text'] for event in events]
    return {'encodings':{name:sum(count(text,name) for text in texts) for name in ENCODINGS},
            'items':len(texts),'basis':'sum of explicit text item encodings; excludes unknown API framing and hidden reasoning',
            'model_mapping_asserted':False,'actual_api_usage':None}
