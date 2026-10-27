import pathlib
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016')
assert (r/'GSE212837_full_validation.json').exists()
assert not (r/'GSM8855042_RDS_download_execution.log').exists()
source=(r/'finish_cosmx_rds.py').read_text()
old="while not (r/'terminal_audit_state.json').exists():time.sleep(20)"
assert old in source
exec(compile(source.replace(old,"assert (r/'GSE212837_full_validation.json').exists()"),str(r/'finish_cosmx_rds.py'),'exec'))
