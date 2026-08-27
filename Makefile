.PHONY: all paper test list rivet check-rivet clean

RIVET_PLUGIN := build/RivetHerwigPolarizedPheno.so
RIVET_SOURCES := \
	analyses/rivet/dis/COMPASS_2009_I820721.cc \
	analyses/rivet/dis/COMPASS_2010_I843494.cc \
	analyses/rivet/dis/COMPASS_2016_I1357198.cc \
	analyses/rivet/dis/COMPASS_2017_I1444985.cc \
	analyses/rivet/dis/COMPASS_2017_I1483098.cc \
	analyses/rivet/dis/COMPASS_2017_I1501480.cc \
	analyses/rivet/dis/HERMES_2007_I726689.cc \
	analyses/rivet/dis/HERMES_2007_I726689_LEGACY.cc \
	analyses/rivet/dis/HERMES_2019_I1698889.cc \
	analyses/rivet/pp/STAR_2019_I1708793.cc \
	analyses/rivet/pp/STAR_2021_I1850855.cc \
	analyses/rivet/pp/STAR_2022_I1949588.cc \
	analyses/rivet/pp/PHENIX_2023_I2033856.cc
RIVET_ANALYSES := \
	COMPASS_2009_I820721 COMPASS_2010_I843494 COMPASS_2016_I1357198 \
	COMPASS_2017_I1444985 COMPASS_2017_I1483098 COMPASS_2017_I1501480 \
	HERMES_2007_I726689 HERMES_2007_I726689_LEGACY HERMES_2019_I1698889 \
	STAR_2019_I1708793 STAR_2021_I1850855 STAR_2022_I1949588 \
	PHENIX_2023_I2033856
RIVET_CXX ?= $(shell \
	configured=$$(rivet-config --cxx 2>/dev/null | awk '{print $$1}'); \
	if test -n "$$configured" && test -x "$$configured"; then printf '%s' "$$configured"; \
	elif command -v g++-16 >/dev/null 2>&1; then command -v g++-16; \
	elif command -v g++-15 >/dev/null 2>&1; then command -v g++-15; \
	elif command -v g++ >/dev/null 2>&1; then command -v g++; \
	else command -v c++; fi)

all: paper

paper: main.pdf

main.pdf: main.tex references.bib jheppub.sty JHEP.bst
	latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

test:
	PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s scripts/tests -p 'test_*.py'

list:
	python3 scripts/run_phenomenology_campaign.py list

rivet: $(RIVET_PLUGIN)

$(RIVET_PLUGIN): $(RIVET_SOURCES) analyses/rivet/dis/COMPASSInclusiveDIS.hh analyses/rivet/dis/COMPASSSIDIS.hh analyses/rivet/dis/COMPASSSIDISBinning.hh analyses/rivet/pp/STARPolarizedJets.hh
	mkdir -p build
	CXX="$(RIVET_CXX)" rivet-build $@ $(RIVET_SOURCES) \
		-I$(CURDIR)/analyses/rivet/dis -I$(CURDIR)/analyses/rivet/pp

check-rivet: rivet
	@RIVET_ANALYSIS_PATH="$(CURDIR)/build:$(CURDIR)/analyses/rivet/dis:$(CURDIR)/analyses/rivet/pp:$${RIVET_ANALYSIS_PATH}"; \
	RIVET_DATA_PATH="$(CURDIR)/analyses/rivet/dis:$(CURDIR)/analyses/rivet/pp:$${RIVET_DATA_PATH}"; \
	export RIVET_ANALYSIS_PATH RIVET_DATA_PATH; \
	for analysis in $(RIVET_ANALYSES); do rivet --show-analysis "$$analysis" >/dev/null || exit 1; done

clean:
	latexmk -C main.tex
