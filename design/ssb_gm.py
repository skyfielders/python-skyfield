from skyfield.data.gravitational_parameters import GM_dict

non_barycenters = [id for id in GM_dict.keys() if not 0 <= id <= 9]
ssb_gm = sum([GM_dict[id] for id in non_barycenters])
print(ssb_gm)
