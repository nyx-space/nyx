from nyx_space.anise import instrument

__all__ = []

for __item__ in dir(instrument):
    if __item__ and __item__[0] != "_":
        print(f"-> adding {__item__}")
        locals()[__item__] = getattr(instrument, __item__)
        __all__ += [__item__]
