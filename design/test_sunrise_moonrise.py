#!/usr/bin/env python
#
# The new 'rising' routine (a) takes only a few steps before returning
# its answer, and (b) ignores the real angular velocity of the target
# across the sky.  It should therefore perform worst in the case of the
# Moon.  Is it able to match USNO predictions for one year?

import argparse
import datetime as dt
import sys
from math import tau
from textwrap import fill
from time import time

from skyfield import almanac
from skyfield.api import load, wgs84
from skyfield.naifcodes import name_codes

# From:
# https://aa.usno.navy.mil/calculated/rstt/year?ID=AA&year=2023&task=0&lat=36.95&lon=-112.52&label=Fredonia%2C+AZ&tz=7&tz_sign=-1&submit=Get+Data

SUN_TABLE = """\
             o  ,    o  ,                                    FREDONIA, AZ                              Astronomical Applications Dept.
Location: W112 31, N36 57                          Rise and Set for the Sun for 2023                   U. S. Naval Observatory
                                                                                                       Washington, DC  20392-5420
                                                      Zone:  7h West of Greenwich


       Jan.       Feb.       Mar.       Apr.       May        June       July       Aug.       Sept.      Oct.       Nov.       Dec.
Day Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set
     h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m
01  0743 1724  0733 1755  0701 1824  0616 1852  0536 1919  0512 1944  0514 1954  0535 1937  0601 1859  0625 1814  0654 1733  0724 1714
02  0743 1725  0732 1756  0700 1825  0615 1853  0535 1920  0512 1945  0515 1953  0536 1936  0602 1857  0626 1812  0655 1732  0725 1713
03  0744 1726  0731 1757  0658 1826  0613 1854  0534 1921  0512 1945  0515 1953  0537 1935  0602 1856  0627 1811  0656 1731  0726 1713
04  0744 1727  0730 1758  0657 1827  0612 1855  0533 1922  0511 1946  0516 1953  0538 1934  0603 1854  0628 1809  0657 1730  0727 1713
05  0744 1728  0729 1759  0656 1828  0610 1856  0532 1922  0511 1946  0516 1953  0539 1933  0604 1853  0629 1808  0658 1729  0728 1713
06  0744 1728  0728 1801  0654 1829  0609 1857  0531 1923  0511 1947  0517 1953  0539 1932  0605 1851  0630 1806  0659 1728  0729 1713
07  0744 1729  0727 1802  0653 1830  0607 1858  0530 1924  0511 1948  0517 1953  0540 1931  0606 1850  0630 1805  0700 1727  0730 1713
08  0744 1730  0726 1803  0651 1831  0606 1859  0529 1925  0510 1948  0518 1952  0541 1930  0607 1848  0631 1803  0701 1726  0731 1713
09  0743 1731  0725 1804  0650 1832  0604 1859  0528 1926  0510 1949  0519 1952  0542 1929  0607 1847  0632 1802  0702 1725  0731 1713
10  0743 1732  0724 1805  0648 1833  0603 1900  0527 1927  0510 1949  0519 1952  0543 1928  0608 1845  0633 1800  0703 1724  0732 1713
11  0743 1733  0723 1806  0647 1834  0602 1901  0526 1928  0510 1950  0520 1951  0544 1926  0609 1844  0634 1759  0704 1724  0733 1714
12  0743 1734  0722 1807  0646 1835  0600 1902  0525 1929  0510 1950  0520 1951  0544 1925  0610 1842  0635 1758  0705 1723  0734 1714
13  0743 1735  0721 1808  0644 1835  0559 1903  0524 1929  0510 1950  0521 1950  0545 1924  0611 1841  0636 1756  0706 1722  0734 1714
14  0743 1736  0720 1809  0643 1836  0557 1904  0523 1930  0510 1951  0522 1950  0546 1923  0611 1839  0637 1755  0707 1721  0735 1714
15  0742 1737  0719 1810  0641 1837  0556 1905  0522 1931  0510 1951  0522 1949  0547 1922  0612 1838  0638 1754  0708 1721  0736 1715
16  0742 1738  0718 1811  0640 1838  0555 1906  0521 1932  0510 1952  0523 1949  0548 1920  0613 1836  0639 1752  0709 1720  0736 1715
17  0742 1739  0716 1812  0638 1839  0553 1906  0521 1933  0510 1952  0524 1948  0549 1919  0614 1835  0639 1751  0710 1719  0737 1715
18  0741 1740  0715 1813  0637 1840  0552 1907  0520 1934  0510 1952  0525 1948  0549 1918  0615 1833  0640 1750  0712 1719  0738 1716
19  0741 1741  0714 1814  0635 1841  0551 1908  0519 1934  0511 1952  0525 1947  0550 1917  0615 1832  0641 1748  0713 1718  0738 1716
20  0740 1742  0713 1815  0634 1842  0549 1909  0518 1935  0511 1953  0526 1947  0551 1915  0616 1830  0642 1747  0714 1717  0739 1716
21  0740 1743  0711 1816  0632 1843  0548 1910  0518 1936  0511 1953  0527 1946  0552 1914  0617 1829  0643 1746  0715 1717  0739 1717
22  0739 1744  0710 1817  0631 1844  0547 1911  0517 1937  0511 1953  0527 1945  0553 1913  0618 1827  0644 1744  0716 1716  0740 1717
23  0739 1745  0709 1818  0629 1844  0546 1912  0517 1938  0511 1953  0528 1945  0554 1911  0619 1826  0645 1743  0717 1716  0740 1718
24  0738 1746  0708 1819  0628 1845  0544 1913  0516 1938  0512 1953  0529 1944  0554 1910  0619 1824  0646 1742  0718 1716  0741 1719
25  0738 1748  0706 1820  0626 1846  0543 1914  0515 1939  0512 1954  0530 1943  0555 1909  0620 1823  0647 1741  0719 1715  0741 1719
26  0737 1749  0705 1821  0625 1847  0542 1914  0515 1940  0512 1954  0531 1942  0556 1907  0621 1821  0648 1740  0720 1715  0742 1720
27  0736 1750  0704 1822  0623 1848  0541 1915  0514 1941  0513 1954  0531 1941  0557 1906  0622 1820  0649 1738  0721 1715  0742 1720
28  0736 1751  0702 1823  0622 1849  0539 1916  0514 1941  0513 1954  0532 1941  0558 1904  0623 1818  0650 1737  0722 1714  0742 1721
29  0735 1752             0620 1850  0538 1917  0513 1942  0513 1954  0533 1940  0558 1903  0624 1817  0651 1736  0723 1714  0743 1722
30  0734 1753             0619 1851  0537 1918  0513 1943  0514 1954  0534 1939  0559 1902  0624 1815  0652 1735  0723 1714  0743 1722
31  0733 1754             0617 1852             0513 1943             0535 1938  0600 1900             0653 1734             0743 1723

Add one hour for daylight time, if and when in use.
"""

MOON_TABLE = """\
             o  ,    o  ,                                    FREDONIA, AZ                              Astronomical Applications Dept.
Location: W112 31, N36 57                         Rise and Set for the Moon for 2023                   U. S. Naval Observatory
                                                                                                       Washington, DC  20392-5420
                                                      Zone:  7h West of Greenwich


       Jan.       Feb.       Mar.       Apr.       May        June       July       Aug.       Sept.      Oct.       Nov.       Dec.
Day Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set
     h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m
01  1342 0257  1414 0456  1259 0345  1442 0425  1529 0345  1726 0323  1839 0306  2011 0508  2015 0749  1940 0901  2025 1113  2111 1128
02  1413 0400  1506 0549  1354 0434  1542 0453  1630 0408  1837 0354  1948 0359  2049 0627  2043 0902  2015 1013  2122 1208  2212 1203
03  1449 0502  1602 0636  1453 0516  1642 0518  1732 0431  1950 0431  2048 0503  2121 0746  2112 1013  2056 1123  2223 1255  2312 1232
04  1530 0603  1701 0716  1553 0552  1742 0542  1837 0456  2101 0518  2138 0617  2149 0900  2144 1123  2144 1228  2324 1333       1257
05  1618 0700  1801 0751  1653 0623  1843 0605  1946 0524  2205 0616  2219 0736  2216 1012  2221 1232  2237 1327       1404  0010 1320
06  1711 0752  1901 0820  1753 0650  1946 0629  2057 0558  2300 0724  2253 0853  2244 1121  2303 1338  2335 1417  0025 1431  0108 1341
07  1809 0837  2000 0846  1852 0715  2052 0655  2208 0638  2345 0839  2322 1007  2313 1230  2352 1439       1459  0124 1455  0206 1403
08  1908 0916  2059 0910  1952 0738  2200 0724  2316 0729       0954  2349 1118  2346 1337       1534  0036 1534  0222 1517  0305 1427
09  2008 0949  2158 0933  2053 0801  2310 0759       0829  0021 1108       1226       1443  0046 1620  0136 1603  0320 1539  0407 1453
10  2107 1017  2259 0956  2155 0825       0842  0015 0938  0052 1219  0015 1333  0023 1547  0145 1659  0236 1629  0419 1602  0512 1524
11  2206 1042       1021  2301 0852  0019 0935  0104 1051  0119 1327  0042 1439  0107 1645  0245 1732  0335 1652  0520 1626  0621 1603
12  2305 1106  0002 1049       0922  0123 1037  0145 1205  0145 1433  0111 1545  0157 1736  0345 1800  0433 1714  0623 1655  0730 1651
13       1129  0108 1122  0009 0959  0218 1147  0219 1317  0211 1539  0145 1650  0253 1820  0445 1825  0531 1736  0730 1728  0838 1749
14  0005 1153  0218 1202  0118 1045  0305 1301  0248 1426  0239 1646  0224 1752  0352 1857  0543 1847  0631 1759  0839 1810  0938 1857
15  0107 1219  0329 1252  0226 1141  0343 1415  0315 1534  0309 1752  0310 1848  0452 1929  0641 1909  0732 1824  0947 1901  1029 2011
16  0213 1249  0437 1354  0328 1246  0416 1527  0341 1642  0345 1857  0402 1938  0552 1956  0739 1931  0836 1854  1050 2001  1111 2126
17  0323 1326  0539 1506  0422 1400  0445 1638  0408 1749  0426 1958  0459 2020  0651 2020  0838 1955  0942 1929  1146 2110  1146 2239
18  0436 1412  0631 1624  0507 1516  0512 1748  0437 1857  0515 2054  0559 2056  0749 2042  0940 2021  1050 2013  1233 2223  1216 2350
19  0549 1509  0714 1743  0544 1632  0539 1857  0509 2004  0609 2141  0659 2126  0847 2104  1044 2052  1156 2106  1311 2336  1243
20  0657 1619  0749 1859  0617 1747  0607 2005  0547 2109  0707 2221  0759 2152  0945 2126  1150 2129  1256 2209  1344       1309 0059
21  0756 1736  0820 2013  0645 1858  0638 2114  0632 2209  0807 2255  0857 2215  1044 2150  1257 2215  1349 2318  1412 0048  1336 0207
22  0844 1856  0848 2123  0713 2009  0713 2221  0722 2301  0907 2323  0955 2237  1146 2218  1402 2312  1433       1439 0158  1405 0316
23  0923 2013  0915 2232  0740 2118  0753 2323  0818 2346  1007 2348  1053 2259  1251 2251  1501       1510 0032  1506 0307  1438 0425
24  0955 2126  0943 2339  0810 2226  0840       0918       1105       1151 2322  1359 2332  1552 0019  1542 0146  1533 0417  1517 0534
25  1024 2236  1013       0842 2333  0933 0020  1018 0023  1203 0012  1252 2348  1508       1635 0133  1610 0300  1604 0528  1603 0641
26  1050 2343  1047 0045  0919       1030 0109  1118 0055  1302 0034  1357       1613 0023  1711 0249  1638 0412  1640 0639  1656 0743
27  1117       1125 0149  1002 0037  1130 0150  1217 0122  1403 0057  1505 0018  1711 0126  1743 0406  1705 0524  1722 0749  1755 0838
28  1144 0048  1209 0249  1050 0137  1230 0225  1316 0146  1506 0121  1616 0055  1801 0239  1812 0522  1735 0636  1812 0856  1857 0923
29  1214 0152             1144 0229  1330 0254  1415 0209  1614 0149  1726 0142  1842 0357  1839 0635  1808 0748  1908 0956  1959 1001
30  1249 0256             1242 0314  1430 0321  1516 0232  1726 0223  1830 0240  1916 0517  1908 0748  1847 0900  2009 1046  2100 1032
31  1328 0357             1342 0352             1619 0256             1926 0350  1946 0634             1933 1010             2159 1058
               Note: Blank spaces in the table indicate that a rising or a setting did not occur during that 24 hr interval.

Add one hour for daylight time, if and when in use.
"""

MOON_TABLE_2 = """
             o  ,    o  ,                                     LATITUDE 70                              Astronomical Applications Dept.
Location:  000 00, N70 00                         Rise and Set for the Moon for 2023                   U. S. Naval Observatory
                                                                                                       Washington, DC  20392-5420
                                                            Universal Time


       Jan.       Feb.       Mar.       Apr.       May        June       July       Aug.       Sept.      Oct.       Nov.       Dec.
Day Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set  Rise  Set
     h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m   h m  h m
01  1047 0408  **** ****  **** ****  **** ****  1345 0421  1917 0100  ---- ----  ---- ----  1939 0619  1625 1026  **** ****  **** ****
02  0950 0641  **** ****  **** ****  1216 0705  1538 0356       0016  ---- ----  2301       1910 0838  **** ****  **** ****  **** ****
03  **** ****  **** ****  **** ****  1425 0626  1734 0332  ---- ----  ---- ----  2214 0416  1837 1057  **** ****  **** ****  1922 1507
04  **** ****  **** ****  **** ****  1620 0558  1943 0306  ---- ----  ---- ----  2143 0656  1743 1337  **** ****  **** ****  2137 1420
05  **** ****  **** ****  1225 0935  1814 0534  2232 0232  ---- ----  ---- ----  2115 0915  **** ****  **** ****  **** ****  2331 1350
06  **** ****  1512 1053  1456 0834  2013 0509       0125  ---- ----  0144 0406  2047 1127  **** ****  **** ****  2208 1635       1325
07  **** ****  1726 1009  1657 0802  2229 0441  ---- ----  ---- ----  0025 0720  2010 1345  **** ****  **** ****       1557  0119 1302
07                                                                    2349
08  **** ****  1921 0941  1850 0736       0401  ---- ----  ---- ----  2322 0942  1853 1646  **** ****  **** ****  0012 1529  0311 1236
09  **** ****  2112 0917  2043 0713  ---- ----  ---- ----  0252 0719  2255 1151  **** ****  **** ****  2221 1903  0204 1505  0513 1202
10  1745 1226  2305 0854  2244 0648  ---- ----  ---- ----  0201 0957  2226 1400  **** ****  **** ****       1806  0354 1440  0752 1059
11  1949 1149       0829       0618  ---- ----  ---- ----  0130 1210  2146 1619  **** ****  **** ****  0046 1733  0550 1411  ---- ----
12  2142 1123  0108 0757  0108 0530  ---- ----  ---- ----  0103 1416  **** ****  **** ****       2022  0244 1707  0803 1331  ---- ----
13  2332 1100  0348 0656  ---- ----  ---- ----  0422 0958  0037 1624  **** ****  **** ****  0108 1940  0435 1643  ---- ----  ---- ----
14       1037  ---- ----  ---- ----  ---- ----  0340 1224  0006 1851  **** ****  **** ****  0318 1911  0626 1617  ---- ----  ---- ----
14                                                         2319
15  0129 1011  ---- ----  ---- ----  0716 0922  0310 1435  **** ****  **** ****       2316  0511 1846  0826 1546  ---- ----  ---- ----
16  0342 0934  ---- ----  ---- ----  0554 1235  0244 1642  **** ****  **** ****  0048 2152  0701 1822  1053 1453  ---- ----  ---- ----
17  ---- ----  ---- ----  ---- ----  0517 1458  0217 1856  **** ****  **** ****  0342 2117  0853 1755  ---- ----  ---- ----  1402 1907
18  ---- ----  ---- ----  ---- ----  0449 1709  0143 2144  **** ****  **** ****  0543 2050  1057 1720  ---- ----  ---- ----  1314 2141
19  ---- ----  1128 1205  0816 1242  0422 1921  0037       **** ****  **** ****  0733 2027  1344 1607  ---- ----  ---- ----  1243 2354
20  ---- ----  0932 1600  0724 1527  0352 2146  **** ****  **** ****  0347 0022  0922 2003  ---- ----  ---- ----  1535 2142  1216
20                                                                         2328
21  ---- ----  0854 1830  0651 1747  0310       **** ****  **** ****  0609 2257  1115 1935  ---- ----  ---- ----  1454       1148 0203
22  ---- ----  0826 2044  0624 1959  **** ****  **** ****  **** ****  0805 2232  1325 1856  ---- ----  ---- ----  1425 0007  1113 0419
23  1149 1628  0759 2254  0557 2216  **** ****  **** ****  0627 0149  0954 2209  ---- ----  ---- ----  1818 2119  1357 0219  1006 0711
24  1059 1910  0731       0523       **** ****  **** ****  0837 0106  1144 2144  ---- ----  ---- ----  1708       1327 0432  **** ****
25  1029 2125  0652 0114  0424 0059  **** ****  **** ****  1030 0038  1341 2114  ---- ----  ---- ----  1632 0020  1246 0657  **** ****
26  1003 2331  **** ****  **** ****  **** ****  0610 0446  1220 0014  1604 2026  ---- ----  1927       1603 0243  **** ****  **** ****
26                                                              2350
27  0938       **** ****  **** ****  **** ****  0906 0320  1413 2324  ---- ----  ---- ----  1838 0033  1534 0459  **** ****  **** ****
28  0908 0138  **** ****  **** ****  **** ****  1108 0244  1620 2250  ---- ----  ---- ----  1806 0315  1500 0719  **** ****  **** ****
29  0822 0402             **** ****  0925 0544  1301 0218  1916 2135  ---- ----  2208       1737 0537  1403 1005  **** ****  **** ****
30  **** ****             **** ****  1146 0452  1454 0154  ---- ----  ---- ----  2046 0029  1706 0757  **** ****  **** ****  1625 1351
31  **** ****             **** ****             1654 0130             ---- ----  2008 0350             **** ****             1901 1244
               Note: Blank spaces in the table indicate that a rising or a setting did not occur during that 24 hr interval.

(**** object continuously above horizon)                                                      (---- object continuously below horizon)

Add one hour for daylight time, if and when in use.
"""

def test(regime, table, e, body, t0, t1, topo, timezone, verbose):
    usno_rises = []

    for line in table.splitlines():
        if not line[0:1].isdigit():
            continue
        day_of_month = int(line[:2])
        months_and_texts = [(i+1, line[4+i*11:8+i*11]) for i in range(0, 12)]
        for month, text in months_and_texts:
            if text.isdigit():
                hours = int(text[0:2])
                minutes = int(text[2:4])
                usno_rises.append((month, day_of_month, hours * 60 + minutes))

    #usno_rises.sort()
    #print(usno_rises)

    observer = e['Earth'] + topo

    if regime == 'old':
        horizon = -0.8333
        if body.target == name_codes['SUN']:
            f = almanac.sunrise_sunset(e, topo)
        else:
            f = almanac.risings_and_settings(e, body, topo, horizon)
            f.step_days = 1 / 1440.0
        tt = time()
        t, y = almanac.find_discrete(t0, t1, f)
        duration = time() - tt
        t = t[y == 1]
    else:
        horizon = -0.8333333333333333
        tt = time()
        t, y = almanac.find_risings(observer, body, t0, t1)
        t = t[y == True]
        duration = time() - tt
    print(f'Compute time: {duration:.3f} seconds')

    # First, before worrying about the USNO: how close did we get to
    # pinning down the actual moment of rising?  Let's compute the
    # target's altitude and see, worse-case, how far it deviates from
    # zero.

    alt, az, distance = observer.at(t).observe(body).apparent().altaz()

    if regime == 'new' and body.target == name_codes['MOON']:
        horizon_array = (almanac._refraction_radians
                         - almanac._moon_radius_m / distance.m)
        horizon_array = horizon_array / tau * 360.0
    else:
        horizon_array = horizon + 0.0 * t.tt

    vs_horizon = alt.degrees - horizon_array

    min_as = min(vs_horizon) * 3600.0
    max_as = max(vs_horizon) * 3600.0

    print(f'Altitude vs horizon: min {min_as:f} arcseconds,',
          f' max {max_as:f} arcseconds')

    def show_shot():
        h = horizon_array[i]
        d = t[i].utc_datetime()
        for offset in -1, 0, 1:
            d2 = d + dt.timedelta(milliseconds=offset)
            t2 = ts.from_datetime(d2)
            alt2, _, _ = observer.at(t2).observe(body).apparent().altaz()
            print('  {}  {}  {: .6f} arcseconds'.format(
                t2.utc_strftime('%Y-%m-%d %H:%M:%S.%f'),
                alt2.dstr(6),
                (alt2.degrees - h) * 3600.0,
            ))
        if verbose:
            def f(t2):
                alt2, _, _ = observer.at(t2).observe(body).apparent().altaz()
                arcseconds = (alt2.degrees - h) * 3600.0
                return arcseconds > 0.0
            hour = 1.0 / 24.0
            #µsec = hour / 3600.00 / 1e6
            f.step_days = hour / 4.0
            t3, y3 = almanac.find_discrete(t[i] - hour, t[i] + hour, f)#, µsec)
            print('  {}  {}  <-- true zero crossing'.format(
                t3[0].utc_strftime('%Y-%m-%d %H:%M:%S.%f'),
                y3,
            ))

    ts = t.ts

    i = vs_horizon.argmin()
    print(f'Worst undershot is at index {i}:')
    show_shot()

    i = vs_horizon.argmax()
    print(f'Worst overshot is at index {i}:')
    show_shot()

    # Okay, now it's time to turn to the USNO table and see how well we
    # stack up again it.

    t += timezone

    skyfield_rises = []

    thirty_seconds = 1.0 / 24.0 / 60.0 / 2.0  # to round to next minute
    for ti in t + thirty_seconds:
        u = ti.utc
        tup = u.month, u.day, u.hour * 60 + u.minute
        skyfield_rises.append(tup)

    usno_dict = {(month, day): time for month, day, time in usno_rises}
    skyfield_dict = {(month, day): time for month, day, time in skyfield_rises}

    errors = []

    for key in sorted(usno_dict.keys() | skyfield_dict.keys()):
        month, day = key
        u = usno_dict.get(key)
        s = skyfield_dict.get(key)
        if not u:
            errors.append('usno-miss')
            print('USNO does not list:', month, day, divmod(int(s), 60))
        elif not s:
            errors.append('skyfield-miss')
            print('Skyfield is missing:', month, day, divmod(u, 60))
        elif u == s:
            errors.append('.')
        else:
            errors.append(f'{u - s}')
    print(fill(''.join(errors), 78))
    print()

def main(argv):
    parser = argparse.ArgumentParser(description='test risings against USNO')
    parser.add_argument('regime', choices=('old', 'new'),
                        help='whether to use old or new rising finder')
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args(argv)
    regime = args.regime
    verbose = args.verbose

    ts = load.timescale()
    e = load('de421.bsp')

    # On 'RuntimeWarning: invalid value encountered in sqrt', stop with
    # a traceback instead of just printing a warning and moving on.
    import warnings
    warnings.filterwarnings("error")

    t0 = ts.utc(2023, 1, 1, 7)
    t1 = ts.utc(2024, 1, 1, 7)
    fredonia = wgs84.latlon(36 + 57/60.0, - (112 + 31/60.0))
    timezone = dt.timedelta(hours=-7)
    test(regime, SUN_TABLE, e, e['Sun'], t0, t1, fredonia, timezone, verbose)
    test(regime, MOON_TABLE, e, e['Moon'], t0, t1, fredonia, timezone, verbose)

    t0 = ts.utc(2023, 1, 1, 0)
    t1 = ts.utc(2024, 1, 1, 0)
    timezone = dt.timedelta(hours=0)
    lat70 = wgs84.latlon(70, 0)
    test(regime, MOON_TABLE_2, e, e['Moon'], t0, t1, lat70, timezone, verbose)

    e.close()

if __name__ == '__main__':
    main(sys.argv[1:])
