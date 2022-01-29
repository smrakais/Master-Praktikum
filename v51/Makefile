all: build/main.pdf

# hier Python-Skripte:
#Z.B
#build/plot1.pdf: plot1.py matplotlibrc header-matplotlib.tex | build
#	TEXINPUTS="$(call translate,$(pwd):)" python plot1.py

# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
#Z.B
#build/main.pdf: build/plot1.pdf

build/main.pdf: FORCE | build
	  TEXINPUTS="$(call translate,build:)" \
	  BIBINPUTS=build: \
	  max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	main.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean