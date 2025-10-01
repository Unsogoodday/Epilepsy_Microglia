# Variables

PYTHON=python
SRC=src
DATA=data
RESULTS=results

# Default target
all: env download preprocess integrate subset analysis figures

# 0. download environments
env:
	$(PYTHON) -m pip install -r requirements.txt

# 1. Download and unpack datasets → data/raw
$(DATA)/raw/done: datasets.txt
	$(PYTHON) $(SRC)/download.py \
		--datasets datasets.txt \
		--outdir $(DATA)/raw
	touch $@

# 2. QC raw data → data/processed
$(DATA)/processed/done: $(DATA)/raw/done
	$(PYTHON) $(SRC)/preprocess.py \
		--input $(DATA)/raw \
		--output $(DATA)/processed
	touch $@

# 3. Integrate
$(DATA)/integrated/done: $(DATA)/processed/done
	$(PYTHON) $(SRC)/integration.py \
		--input $(DATA)/processed \
		--output $(DATA)/integrated
	touch $@

# 4. Subset
$(DATA)/subset/done: $(DATA)/integrated/done
	$(PYTHON) $(SRC)/subset.py \
		--input $(DATA)/integrated \
		--output $(DATA)/subset
	touch $@

# 5. Analysis (integrated)
$(RESULTS)/tables/integrated_done: $(DATA)/integrated/done
	$(PYTHON) $(SRC)/analysis_integrated.py \
		--input $(DATA)/integrated \
		--output $(RESULTS)/tables
	touch $@

# 6. Analysis (subset)
$(RESULTS)/tables/subset_done: $(DATA)/subset/done
	$(PYTHON) $(SRC)/analysis_subset.py \
		--input $(DATA)/subset \
		--output $(RESULTS)/tables
	touch $@

# 7. Figures
$(RESULTS)/figures/done: $(RESULTS)/tables/integrated_done $(RESULTS)/tables/subset_done
	$(PYTHON) $(SRC)/plotting.py \
		--input $(RESULTS)/tables \
		--output $(RESULTS)/figures
	touch $@


# Aliases
preprocess: $(DATA)/processed/done
integrate: $(DATA)/integrated/done
subset: $(DATA)/subset/done
analysis: $(RESULTS)/tables/integrated_done $(RESULTS)/tables/subset_done
figures: $(RESULTS)/figures/done

# Utility
clean:
	rm -rf $(DATA)/* $(RESULTS)/*

