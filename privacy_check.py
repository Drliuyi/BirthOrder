#!/usr/bin/env python3
"""Fail-closed code-only check of working, staged or committed files.

This is a release safeguard, not proof that arbitrary future edits are anonymous.
Manual review is still required before adding a new file to release_files.json.
"""
from pathlib import Path
import argparse, ast, csv, hashlib, io, json, re, subprocess

ROOT=Path(__file__).resolve().parent
MANIFEST='reference/release_files.json'
DICTIONARIES={
    'reference/diseases.csv': {'phecode','description','group','groupnum','color'},
    'reference/sex_restrictions.csv': {'phecode','male_only','female_only'},
    'reference/mediators.csv': {'indicator_id','category','original_label','type','measurement_instance'},
    'reference/table1.csv': {'column','label','type'},
}
FORBIDDEN_PATH=re.compile(r'(^|/)(data|input|inputs|results|output|logs|cache|checkpoints|private|draws|\.ssh|\.env)(/|$)',re.I)
CONTENT_PATTERNS={
    'personal machine path': re.compile(r'/(?:Users|home|Volumes)/[A-Za-z0-9_.-]+/'),
    'private network address': re.compile(r'\b(?:10(?:\.\d{1,3}){3}|192\.168(?:\.\d{1,3}){2}|172\.(?:1[6-9]|2\d|3[01])(?:\.\d{1,3}){2})\b'),
    'GitHub token': re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,})\b'),
    'private key': re.compile(r'-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----'),
    'literal participant identifier': re.compile(r'''(?:\beid\b|participant_id|sample_eid)\s*[:=]\s*["']?\d{6,8}\b''',re.I),
    'embedded participant rows': re.compile(r'''["'](?:eid|participant_id)["']\s*:\s*\[\s*["']?\d{6,8}''',re.I),
}

def check_blob(name,blob,expected):
    errors=[]
    path=Path(name)
    if name not in expected:return ['file not explicitly reviewed/allowlisted']
    if FORBIDDEN_PATH.search(name):errors.append('forbidden data/cache directory')
    if b'\0' in blob:return errors+['binary content']
    try:text=blob.decode('utf-8')
    except UnicodeDecodeError:return errors+['not UTF-8 source']
    if len(blob)>3_000_000:errors.append('unexpectedly large source file')
    for label,pattern in CONTENT_PATTERNS.items():
        if pattern.search(text):errors.append(label)
    if path.suffix=='.csv':
        if name not in DICTIONARIES:errors.append('CSV is not an approved variable dictionary')
        else:
            reader=csv.DictReader(io.StringIO(text))
            fields=set(reader.fieldnames or [])
            if not fields <= DICTIONARIES[name]:errors.append('unexpected dictionary fields')
            for row in reader:
                if any(re.fullmatch(r'\d{6,8}',str(v or '')) for v in row.values()):
                    errors.append('possible individual-level identifier in dictionary');break
    if path.suffix=='.py':
        try:ast.parse(text)
        except SyntaxError:errors.append('invalid Python syntax')
    actual=hashlib.sha256(blob).hexdigest()
    if name!=MANIFEST and actual!=expected[name]['sha256']:errors.append('content changed since reviewed release manifest')
    return errors

def git(*args):
    return subprocess.check_output(['git',*args],cwd=ROOT)

def main():
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--staged',action='store_true')
    group.add_argument('--tracked',action='store_true')
    opts=parser.parse_args()
    manifest=json.loads((ROOT/MANIFEST).read_text())
    expected={r['file']:r for r in manifest['files']};expected[MANIFEST]={}
    if opts.staged:names=git('ls-files','--cached','-z').decode().split('\0')
    elif opts.tracked:names=git('ls-tree','-r','--name-only','-z','HEAD').decode().split('\0')
    else:names=[str(p.relative_to(ROOT))for p in ROOT.rglob('*')if p.is_file()and '.git'not in p.relative_to(ROOT).parts]
    errors={}
    for name in filter(None,names):
        p=ROOT/name
        if p.is_symlink():errors[name]=['symbolic link'];continue
        if opts.staged:blob=git('show',':'+name)
        elif opts.tracked:blob=git('show','HEAD:'+name)
        else:blob=p.read_bytes()
        problem=check_blob(name,blob,expected)
        if problem:errors[name]=problem
    present=set(filter(None,names))
    for name in expected.keys()-present:errors[name]=['reviewed release file missing']
    print(json.dumps({'files_checked':len(present),'violations':errors,'code_only':not errors},indent=2))
    raise SystemExit(bool(errors))

if __name__=='__main__':main()
