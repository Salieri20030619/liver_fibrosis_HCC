"""Append exactly one validated chunk and verify persisted bytes before accepting it."""
import hashlib,os
def append_verified(body, destination, offset, expected_length):
    actual=body.stat().st_size
    if actual!=expected_length:raise ValueError('Chunk length changed before append')
    if (destination.stat().st_size if destination.exists() else 0)!=offset:
        raise ValueError('Destination offset mismatch; preserve evidence')
    h=hashlib.sha256()
    with body.open('rb') as src,destination.open('ab') as dst:
        while True:
            data=src.read(1024*1024)
            if not data:break
            h.update(data);dst.write(data)
        dst.flush();os.fsync(dst.fileno())
    return verify_readback(destination,offset,expected_length,h.hexdigest())
def verify_readback(destination,offset,length,expected_sha256):
    h=hashlib.sha256();remaining=length
    with destination.open('rb') as f:
        f.seek(offset)
        while remaining:
            data=f.read(min(1024*1024,remaining))
            if not data:raise ValueError('Truncated persisted chunk')
            h.update(data);remaining-=len(data)
    observed=h.hexdigest()
    if observed!=expected_sha256:raise ValueError('Persisted chunk SHA256 mismatch; preserve body and destination')
    return {'offset':offset,'bytes':length,'body_sha256':expected_sha256,'readback_sha256':observed,'readback_matches':True}
