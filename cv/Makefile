CV=AlexGomez-CV
OUT=../assets/$(CV).pdf

all: $(OUT)

$(OUT): $(CV).tex
	mkdir -p build
	pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build $(CV).tex
	pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build $(CV).tex
	cp build/$(CV).pdf $(OUT)

clean:
	rm -rf build
