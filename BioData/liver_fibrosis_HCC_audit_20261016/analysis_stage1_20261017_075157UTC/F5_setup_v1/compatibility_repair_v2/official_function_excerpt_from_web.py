def _bind_property(cls, pd_cls, attr, min_version=None):
    def func(self):
        return self._property_map(attr)
    func.__name__ = attr
    func.__qualname__ = f"{cls.__name__}.{attr}"
    original_prop = getattr(pd_cls, attr)
    if isinstance(original_prop, property):
        method = original_prop.fget
    elif isinstance(original_prop, functools.cached_property):
        method = original_prop.func
    else:
        method = original_prop
    try:
        func.__wrapped__ = method
    except Exception:
        pass
    setattr(cls, attr, property(derived_from(pd_cls, version=min_version)(func)))

