# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build
IMAGEDIRS     = GUI Validations

# SVG to PDF conversion
SVG2PDF       = inkscape
SVG2PDF_FLAGS =

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Pattern rule for converting SVG to PDF
%.pdf : %.svg
    $(SVG2PDF) -f $< -A $@

# Build a list of SVG files to convert to PDFs
PDFs := $(foreach dir, $(IMAGEDIRS), $(patsubst %.svg,%.pdf,$(wildcard $(SOURCEDIR)/$(dir)/*.svg)))

# Make a rule to build the PDFs
images: $(PDFs)

clean:
    -rm -rf $(BUILDDIR)/*
    -rm $(PDFs)

latex: images
    # [...]
latexpdf: images
    # [...]

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

