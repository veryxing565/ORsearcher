from __future__ import annotations
import hashlib,json,logging,os,signal,subprocess,time
from pathlib import Path

LOG=logging.getLogger('orsearcher2')
PROCESS_GROUP_SUPPORTED=os.name=='posix'

def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,default=str).encode()).hexdigest()

def file_digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def metadata(path):
    p=Path(path).resolve();s=p.stat()
    return dict(path=str(p),size=s.st_size,mtime_ns=s.st_mtime_ns)

def save(path,value,compact=False):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(value,indent=None if compact else 2,
        separators=(',',':') if compact else None,default=str)+'\n',encoding='utf-8');os.replace(temp,path)

def execute(cmd,output,timeout=None,inputs=(),resume=True,kill_process_group=False):
    """Cache only successful commands; failed/partial output is never a hit."""
    LOG.debug('Command: %s | timeout=%s | output=%s',list(map(str,cmd)),timeout,output)
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    marker=Path(str(output)+'.status.json')
    key=digest(dict(command=list(map(str,cmd)),inputs=[metadata(p) for p in inputs]))
    if resume and output.exists() and marker.exists():
        old=json.loads(marker.read_text())
        if old.get('fingerprint')==key and old.get('status')=='completed':return old
    start=time.monotonic();temp=Path(str(output)+'.partial');status='completed';error='';code=None
    with temp.open('w',encoding='utf-8') as out,Path(str(output)+'.stderr.log').open('w',encoding='utf-8') as err:
        try:
            if kill_process_group and PROCESS_GROUP_SUPPORTED:
                # MetaEuk starts shell/MMseqs descendants. Killing only its
                # launcher on timeout would leave CPU-consuming orphans.
                with subprocess.Popen(list(map(str,cmd)),stdout=out,stderr=err,start_new_session=True) as p:
                    try:code=p.wait(timeout=timeout)
                    except subprocess.TimeoutExpired:
                        try:os.killpg(p.pid,signal.SIGTERM)
                        except ProcessLookupError:pass
                        try:p.wait(timeout=2)
                        except subprocess.TimeoutExpired:pass
                        # Also stop any descendant ignoring TERM, even if the
                        # launcher already exited. This is our own session only.
                        try:os.killpg(p.pid,signal.SIGKILL)
                        except ProcessLookupError:pass
                        p.wait();raise
            else:
                p=subprocess.run(list(map(str,cmd)),stdout=out,stderr=err,timeout=timeout)
                code=p.returncode
            if code:status='error';error=f'exit {code}'
        except subprocess.TimeoutExpired:status='timeout';error=f'limit {timeout}s'
    record=dict(command=list(map(str,cmd)),fingerprint=key,status=status,error=error,exit_code=code,seconds=round(time.monotonic()-start,2),process_group_cleanup=bool(kill_process_group and PROCESS_GROUP_SUPPORTED))
    if status=='completed':os.replace(temp,output)
    save(marker,record)
    if status!='completed':raise RuntimeError(f'{status}: {cmd[0]}: {error}: {marker}')
    return record

def fasta_if_any(path,reader):
    return list(reader(path)) if Path(path).exists() and Path(path).stat().st_size else []
