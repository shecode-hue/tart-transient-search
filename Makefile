.PHONY: help install patches test run download build search clean

CONFIG ?= config/za-hammanskraal.yaml

help:
	@echo "make install    install dependencies"
	@echo "make patches    apply the upstream tart2ms/tartcargo patches"
	@echo "make test       run the convention regression tests"
	@echo "make run        full pipeline        (CONFIG=$(CONFIG))"
	@echo "make download   fetch HDF only"
	@echo "make build      measurement set + calibration only"
	@echo "make search     fit, peel, image, search only"
	@echo "make clean      remove runs/"

install:
	pip install -r requirements.txt
	pip install -e .

patches:
	python3 scripts/apply_patches.py

test:
	PYTHONPATH=src pytest tests/ -v

run:
	PYTHONPATH=src python3 -m tart_transient run --config $(CONFIG)

download:
	PYTHONPATH=src python3 -m tart_transient download --config $(CONFIG)

build:
	PYTHONPATH=src python3 -m tart_transient build --config $(CONFIG)

search:
	PYTHONPATH=src python3 -m tart_transient search --config $(CONFIG)

clean:
	rm -rf runs/
