import spiceypy as sp

sp.furnsh('de438s.bsp')
sp.furnsh('naif0012.tls')
sp.furnsh('pck00008.tpc')

et = sp.utc2et('2022-11-30 12:00:00')

moon_iau_earth = sp.spkpos('moon', et, 'iau_earth', 'none', 'earth')[0]
moon_eme2k = sp.spkpos('moon', et, 'j2000', 'none', 'earth')[0]

print(f'moon_iau_earth = {moon_iau_earth}\nmoon_eme2k = {moon_eme2k}')


def sxform_val(inf, outf, et):
    print(f"\n{inf} -> {outf}")
    print(sp.sxform(inf, outf, et))


sxform_val('iau_earth', 'j2000', et)
sxform_val('j2000', 'iau_earth', et)

sxform_val('iau_earth', 'iau_mars', et)
sxform_val('iau_mars', 'iau_earth', et)
