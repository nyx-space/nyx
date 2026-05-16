from nyx_space.anise import instrument

__all__ = []

for el in dir(instrument):
    if el and el[0] != "_":
        print(f"-> adding {el}")
        locals()[el] = getattr(instrument, el)
        __all__ += [el]
