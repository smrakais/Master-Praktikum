all: build/protokoll.pdf


build/protokoll.pdf: protokoll.tex data		
	latexmk --lualatex protokoll.tex

data:# build/EMPTY.tex

#ERZEUGE				ZUTATEN		REZEPT(TERMINAL BEFEHL)
build/protokoll.pdf: #build/wuerfel1.tex build/wuerfel4_erg.tex build/cäsium.pdf
build/protokoll.pdf: build
          max_print_line=1048576 \
		  BIBINPUTS=build: \
        latexmk \
          --lualatex \
          --output-directory=build \
          --interaction=nonstopmode \
          --halt-on-error \
        protokoll.tex
#
#build/wuerfel1.tex: make_tables.py | build
#	python make_tables.py
#
#build/wuerfel1_erg.tex: aluminium_i0.py | build
#	python aluminium_i0.py
#
#build/wuerfel4_erg.tex: leastsquares.py | build
#	python leastsquares.py
#
#build/cäsium.pdf: caesium_spektrum.py | build
#	python caesium_spektrum.py

build: 
	mkdir -p build

clean:
	rm -r build
