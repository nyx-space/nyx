from nyx_space.anise import analysis

__all__ = []

for el in dir(analysis):
    if el and el[0] != "_":
        locals()[el] = getattr(analysis, el)
        __all__ += [el]
