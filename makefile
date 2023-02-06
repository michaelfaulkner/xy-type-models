all: multivalued-electrolyte elementary-electrolyte hxy-ecmc hxy-uniform-noise-metropolis \
		hxy-gaussian-noise-metropolis xy-ecmc xy-uniform-noise-metropolis xy-gaussian-noise-metropolis 3dxy

elementary-electrolyte:
	cd src/electrolyte/elementary && $(MAKE)

multivalued-electrolyte:
	cd src/electrolyte/multivalued && $(MAKE)

hxy-ecmc:
	cd src/xy_models/hxy/ecmc && $(MAKE)

hxy-uniform-noise-metropolis:
	cd src/xy_models/hxy/metropolis/uniform_noise_metropolis && $(MAKE)

hxy-gaussian-noise-metropolis:
	cd src/xy_models/hxy/metropolis/gaussian_noise_metropolis && $(MAKE)

xy-ecmc:
	cd src/xy_models/xy/ecmc && $(MAKE)

xy-uniform-noise-metropolis:
	cd src/xy_models/xy/metropolis/uniform_noise_metropolis && $(MAKE)

xy-gaussian-noise-metropolis:
	cd src/xy_models/xy/metropolis/gaussian_noise_metropolis && $(MAKE)

3dxy:
	cd src/3dxy && $(MAKE)
