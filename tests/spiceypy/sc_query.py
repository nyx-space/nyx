import spiceypy as sp

sp.furnsh('de438s.bsp')
sp.furnsh('naif0012.tls')
sp.furnsh('pck00008.tpc')
sp.furnsh('mysc-Moon-J2k.bsp')
# sp.furnsh('mysc-EME2k.bsp')


def print_state(epoch):
    et = sp.utc2et(epoch)
    print()
    print(epoch, et)
    sc_from_moon = sp.spkez(-10003001, et, 'j2000', 'none', 301)[0]
    print("Moon J2000", sc_from_moon)
    sc_from_earth = sp.spkez(-10003001, et, 'j2000', 'none', 399)[0]
    print("EME2000", sc_from_earth)
    moon_from_earth = sp.spkez(301, et, 'j2000', 'none', 399)[0]
    print("Moon from Earth:", moon_from_earth)
    moon_from_emb = sp.spkez(301, et, 'j2000', 'none', 3)[0]
    print("Moon from EMB:", moon_from_emb)
    earth_from_emb = sp.spkez(399, et, 'j2000', 'none', 3)[0]
    print("Earth from EMB:", earth_from_emb)
    sum = -earth_from_emb + moon_from_emb + sc_from_moon
    print("Sum: -earth_from_emb + moon_from_emb + sc_from_moon\n\t", sum)
    print("Delta:", sum - sc_from_earth)


almost_epoch = '29 Nov 2022 06:47:05.000'
# GMAT state: 29 Nov 2022 06:47:28.000,-274550.832878,197503.896653,120515.732759,7.207074,-2.499496,-2.009755
print_state(almost_epoch)

almost_epoch = '2022-12-04 11:59:51.884'
# GMAT state: 04 Dec 2022 12:00:28.884,482.6690997,-1202.72242593,1766.54770503,0.575668890281,-0.97502959975,-1.95884508326
print_state(almost_epoch)
