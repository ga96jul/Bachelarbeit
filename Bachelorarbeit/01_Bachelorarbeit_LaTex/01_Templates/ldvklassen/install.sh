#!/bin/sh

SRCPATH=$(dirname $0)

# Ziel ist das Verzeichnis für lokale TeX-Pakete
TEXMFLOCAL=$(kpsexpand '$TEXMFLOCAL')
LYXDIR=/usr/share/lyx/layouts

# Sind die zu installierenden Dateien vorhanden?
# Teste einzelne Dateien in jedem Unterverzeichnis.
if [ ! -f "$SRCPATH/tex/latex/ldv/ldvbook.cls" -o ! -f "$SRCPATH/doc/latex/ldv/ldvguide.pdf" -o ! -f "$SRCPATH/src/latex/ldv/ldvcommon.dtx" -o ! -f "$SRCPATH/lyx/ldvbook.layout" ]
then
		echo "Konnte die zu installierenden Dateien nicht finden."
		exit 1
fi

# Falls das Verzeichnis noch nicht existiert, lege es an.
if [ ! -d "$TEXMFLOCAL" ]
then
		mkdir -p "$TEXMFLOCAL"
fi

# Die LaTeX-Dateien kopieren. Die Verzeichnisstruktur passt bereits.
echo "Installiere die LaTeX-Dateien"
mkdir -p "$TEXMFLOCAL/tex/latex/ldv"
cp "$SRCPATH/tex/latex/ldv"/*.cls "$SRCPATH/tex/latex/ldv"/*.pdf "$SRCPATH/tex/latex/ldv"/*.eps "$TEXMFLOCAL/tex/latex/ldv"
mkdir -p "$TEXMFLOCAL/source/latex/ldv"
cp "$SRCPATH/src/latex/ldv"/* "$TEXMFLOCAL/source/latex/ldv"
mkdir -p "$TEXMFLOCAL/doc/latex/ldv"
cp "$SRCPATH/doc/latex/ldv"/* "$TEXMFLOCAL/doc/latex/ldv"
mkdir -p "$TEXMFLOCAL/bibtex/bst/ldv/"
cp "$SRCPATH/bibtex/bst/ldv"/*.bst "$TEXMFLOCAL/bibtex/bst/ldv/"

# Die neuen Dateien der TeX-Datenbank bekannt machen
texhash

# Die LyX-Dateien kopieren.
if [ -d "$LYXDIR" ]
then
		echo "Installiere die LyX-Dateien"
		cp -R "$SRCPATH/lyx"/* "$LYXDIR"
fi
