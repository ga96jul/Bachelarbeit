% 2008 LNT Thesis class
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{LNTthesis}[2009/08/24 LNT Thesis class]

    \DeclareOption{DIV12}{\PassOptionsToClass{\CurrentOption}{scrreprt}}
    \DeclareOption{BCOR12mm}{\PassOptionsToClass{\CurrentOption}{scrreprt}}
    \DeclareOption{headsepline}{\PassOptionsToClass{\CurrentOption}{scrreprt}}
    \DeclareOption{twoside}{\PassOptionsToClass{\CurrentOption}{scrreprt}}
    \DeclareOption{11pt}{\PassOptionsToClass{\CurrentOption}{scrreprt}}
    \DeclareOption{a4paper}{\PassOptionsToClass{\CurrentOption}{scrreprt}}
    \DeclareOption{bibtotoc}{\PassOptionsToClass{\CurrentOption}{scrreprt}}
	\DeclareOption{openright}{\PassOptionsToClass{\CurrentOption}{scrreprt}}
	\DeclareOption{cleardoubleempty}{\PassOptionsToClass{\CurrentOption}{scrreprt}}
% If the document passes options to this class, pass them to the scrreprt class
    \DeclareOption*{\PassOptionsToClass{\CurrentOption}{scrreprt}}
    \ExecuteOptions{a4paper,DIV12,BCOR12mm,headsepline,twoside,11pt,bibtotoc,openright,cleardoubleempty}
% Process options (MUST BE BEFORE LOADCLASS COMMAND!)
    \ProcessOptions\relax

% This class is built upon the KOMAscript scrreprt class
    \LoadClass{scrreprt}
    % Parameters of the scrreprt class:
    % DIV12       = document layout: DIV factor = 12 (larger number creates larger pages)
    % BCOR12mm    = binding correction: 12mm
    % headsepline = separate page head by a line
    % twoside     = twosided document
    % 11pt        = font size 11 point
    % openright   = start new chapters only on right pages (uneven numbered pages)
    % more information: http://tug.ctan.org/tex-archive/macros/latex/contrib/koma-script/scrguien.pdf


% Load necessary packages
    \RequirePackage[english]{babel} % english language / Englische Sprache
    \RequirePackage{graphicx}       % for using eps images / Einbinden von EPS Grafiken
    \RequirePackage{verbatim}       % for quickly commenting out large parts of your text / Um viel Text schnell auskommentieren zu koennen
    \RequirePackage{amssymb}        % additional math symbols / Zusaetzliche mathematische Symbole
    \RequirePackage{amsmath}        % additional math commands / Zusaetzliche mathematische Befehle
    \RequirePackage{amsxtra}        % even more math symbols / Noch mehr mathematische Symbole
    \RequirePackage{amsthm}         % theorem environment etc / Theorem Umgebung usw
    \RequirePackage{setspace}		% set line spacing globally / Zeilenabstand global setzen 
    % more information on amsmath: http://www.ctan.org/get/macros/latex/required/amslatex/math/amsldoc.pdf
    \RequirePackage{psfrag}         % psfrag: http://www.ctan.org/get/macros/latex/contrib/psfrag/pfgguide.pdf
    \RequirePackage{subfigure}      % enable subfigures / Ermoeglicht Subfigures (mehrere Figures neben/untereinander)

% Set chapter, section and subsection to regular LaTeX font
    \addtokomafont{sectioning}{\rmfamily}
    \addtokomafont{descriptionlabel}{\rmfamily}

% Set general pagesyle / Allgemeinen Pagestyle setzen
    \pagestyle{headings}

% Revome paragraph indent / Absatzeinzug entfernen
    \parindent 0em
    
% Set line spacing globally / Zeilenabstand global setzen
 	\onehalfspacing

% Set depth of table of contents to 2 (this includes sections and subsections)
 	\setcounter{tocdepth}{2}

% Title page / Titelseite
\newcommand{\LNTtitle}[6]{
	\begingroup
	\onehalfspacing
    \titlehead{
        \vspace*{-40pt}
        \begin{tabular}{llr}
            \includegraphics[height=16.5mm]{TUMlogo_tumcd} &
            \parbox[b]{85mm}{
                    {\large \textsc{Technische Universit\"at M\"unchen} }\\
                    {\large \textsc{Lehrstuhl f\"ur Nachrichtentechnik} }\\
                    {\normalsize Prof.~Dr.~sc.~techn.~Gerhard~Kramer\vspace*{1mm}}
            } &
            \includegraphics[height=17mm]{LNTlogo_tumcd} \\[0.5cm]
            \hline
        \end{tabular}
    }
    \subject{\vspace*{2.5cm}#1}
    \title{#2}
    \author{\vspace*{1cm}\\Vorgelegt von:\\#3}
    \publishers{\vspace*{2.5cm}Betreut von:\\#4\\[-2cm]}
    \uppertitleback{#1 am\\
                Lehrstuhl f\"ur Nachrichtentechnik (LNT)\\
                der Technischen Universit\"at M\"unchen (TUM)\\
                Titel : #2\\
                Autor : #3}
    \lowertitleback{#6}
    %\dedication{#?}
    \date{#5}
    \maketitle
    \endgroup
}

% *******************************************************
% Rechtliche Erkl\"arung

\newcommand{\LNTrecht}[2]{
	\begingroup
	\onehalfspacing
    \thispagestyle{empty}
    \
    \vfill
    Ich versichere hiermit wahrheitsgem\"a\ss, die Arbeit bis auf die dem
    Aufgabensteller bereits bekannte Hilfe selbst\"andig angefertigt, alle
    benutzten Hilfsmittel vollst\"andig und genau angegeben und alles kenntlich
    gemacht zu haben, was aus Arbeiten anderer unver\"andert oder mit
    Ab\"anderung entnommen wurde.\\
    \begin{center}
    \raggedright #2\\
    \vspace*{-2ex}
    \dotfill\\
    Ort, Datum \hfill (#1)\\
    \end{center}
    \endgroup
}   