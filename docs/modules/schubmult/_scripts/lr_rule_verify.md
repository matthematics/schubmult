<a id="schubmult._scripts.lr_rule_verify"></a>

# schubmult.\_scripts.lr\_rule\_verify

<a id="schubmult._scripts.lr_rule_verify.json_key_rep"></a>

#### json\_key\_rep

```python
def json_key_rep(k)
```

Serialize an arbitrary Python object `k` to a JSON-key-safe string.
Uses pickle + base64 with a 'PICKLE:' prefix; falls back to repr.

<a id="schubmult._scripts.lr_rule_verify.json_key_load"></a>

#### json\_key\_load

```python
def json_key_load(s)
```

Inverse of json_key_rep.

<a id="schubmult._scripts.lr_rule_verify.safe_save_recording"></a>

#### safe\_save\_recording

```python
def safe_save_recording(obj, filename, meta=None)
```

Save {meta: {...}, records: {json_key: value}} atomically.

<a id="schubmult._scripts.lr_rule_verify.safe_load_recording"></a>

#### safe\_load\_recording

```python
def safe_load_recording(filename)
```

Load the verification JSON.
Returns (records_dict, meta_dict).
Backwards compatible with old flat mapping format (treats file as records and meta={}).

<a id="schubmult._scripts.lr_rule_verify.recording_saver"></a>

#### recording\_saver

```python
def recording_saver(shared_recording_dict,
                    lock,
                    verification_filename,
                    stop_event,
                    meta=None,
                    sleep_time=10)
```

Background process that periodically snapshots `shared_recording_dict` and
writes it to `verification_filename`. The lock is held only briefly to
snapshot keys; the potentially slow file I/O and proxy lookups are done
outside the lock so workers are not blocked.

<a id="schubmult._scripts.lr_rule_verify.queue_producer"></a>

#### queue\_producer

```python
def queue_producer(task_queue, perms, n, num_processors, skip_id, irreducible,
                   base_level, sep_descs, dominant_only, w0_only,
                   shared_recording_dict, molevsagan)
```

Producer with periodic garbage collection.

