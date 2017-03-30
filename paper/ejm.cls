%%
%% This is file `ejm.cls'
%%
%% CUP European Journal of Applied Mathematics document class
%% Copyright 2001 Cambridge University Press
%%
%% by Mark A. Reed
%%   based on EJM.STY v0.4.
%% updated by David Tranah 29 June 2010
%%
%% Incorporating parts of authordate.sty, by David Rhead.
%% Incorporating subeqn.sty, by Stephen Gildea.
%% Incorporating subeqnarray.sty, by Johannes Braams.
%%
%% Bugs (in the case of unchanged files) should be reported to
%% texline@cup.cam.ac.uk
%%
%% This software may only be used in the preparation of journal articles
%% or books or parts of books to be published by Cambridge University Press.
%% Any other use constitutes an infringement of copyright.
%%
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
%%

\NeedsTeXFormat{LaTeX2e}[1997/12/01]
\ProvidesClass{ejm}[2010/06/29 v1.04 European Journal of Applied Mathematics
  ^^Jdocument class]

\newif\ifprodtf
\newif\ifreferee

\DeclareOption{oneside}{\relax}
\DeclareOption{twoside}{\@twosidetrue \@mparswitchtrue}
\DeclareOption{draft}{\setlength\overfullrule{5\p@}}
\DeclareOption{final}{\setlength\overfullrule{\z@}}
\DeclareOption{onecolumn}{\@twocolumnfalse}
\DeclareOption{twocolumn}{\relax}
\DeclareOption{titlepage}{\relax}
\DeclareOption{notitlepage}{\relax}
\DeclareOption{leqno}{\relax}
\DeclareOption{fleqn}{\relax}
\DeclareOption{prodtf}{\prodtftrue}
\DeclareOption{referee}{\refereetrue}

\ExecuteOptions{twoside,final,onecolumn}
\ProcessOptions\relax

\newif\ifCUPmtlplainloaded
\ifprodtf
  \CUPmtlplainloadedtrue
  \RequirePackage{CUPTimes,ejm2esym}
\fi

\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}
\renewcommand\baselinestretch{}

\newcommand\tusp{\kern1.25\p@\relax}
\newcommand\fusp{\kern2.75\p@\relax}

\renewcommand\normalsize{%
  \@setfontsize\normalsize\@xpt{13}%
  \abovedisplayskip 6.5\p@ \@plus 1\p@ \@minus 1\p@
  \belowdisplayskip \abovedisplayskip
  \abovedisplayshortskip 3\p@ \@plus 1\p@
  \belowdisplayshortskip \abovedisplayshortskip
}

\normalsize

\AtBeginDocument{%
  \ifuseAMSsubequ \else
    \global\let\subequations\CUPsubequations
    \global\let\endsubequations\endCUPsubequations
  \fi
}

\newcommand\small{%
  \@setfontsize\small\@ixpt\@xpt
  \abovedisplayskip 6\p@ \@plus 1\p@ \@minus 1\p@
  \belowdisplayskip \abovedisplayskip
  \abovedisplayshortskip 3\p@ \@plus 1\p@
  \belowdisplayshortskip \abovedisplayshortskip
}

\newcommand\footnotesize{\@setfontsize\footnotesize\@ixpt\@xipt}
\newcommand\scriptsize{\@setfontsize\scriptsize\@viipt\@viiipt}
\newcommand\tiny{\@setfontsize\tiny\@vpt\@vipt}
\newcommand\large{\@setfontsize\large\@xipt{13}}
\newcommand\Large{\@setfontsize\Large\@xivpt{18}}
\ifprodtf
  \newcommand\LARGE{\@setfontsize\LARGE{18}{21}}%
\else
  \newcommand\LARGE{\@setfontsize\LARGE\@xviipt{21}}%
\fi
\newcommand\huge{\@setfontsize\huge\@xxpt{25}}
\newcommand\Huge{\@setfontsize\Huge\@xxvpt{30}}

\newcommand\abstractsize{\@setfontsize\abstractsize\@ixpt\@xiipt}
\newcommand\authorsize{\@setfontsize\authorsize\@viiipt\@xpt}
\newcommand\catchlinesize{\@setfontsize\catchlinesize\@viiipt\@ixpt}
\newcommand\figcapsize{\footnotesize}
\newcommand\bibliosize{\footnotesize}
\newcommand\tablebdysize{\footnotesize}
\newcommand\quotesize{\footnotesize}

\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}

\edef\r@{\ifprodtf mtr\else cmr\fi}
\edef\s@{\ifprodtf mtss\else cmss\fi}

\newcommand\qsmall{%
  \usefont{OT1}{\r@}{m}{n}\fontsize{9.75}{\f@baselineskip}\selectfont}
\newcommand\qsmit{%
  \usefont{OT1}{\r@}{m}{it}\fontsize{9.75}{\f@baselineskip}\selectfont}
\newcommand\qsbf{%
  \usefont{OT1}{\r@}{bx}{n}\fontsize{9.75}{\f@baselineskip}\selectfont}
\newcommand\qssf{%
  \usefont{OT1}{\s@}{m}{n}\fontsize{9.75}{\f@baselineskip}\selectfont}
\newcommand\boldit{%
  \usefont{OT1}{\r@}{bx}{it}\fontsize{10}{\f@baselineskip}\selectfont}

\newcommand\qit{\qsmit \fontsize{\f@size}{10.75}\selectfont}
\newcommand\qbf{\qsbf  \fontsize{\f@size}{10.75}\selectfont}
\newcommand\qsf{\qssf  \fontsize{\f@size}{10.75}\selectfont}
\newcommand\qrm{\qsmall\fontsize{\f@size}{10.75}\selectfont}

\ifprodtf\else
  \DeclareTextFontCommand\textsfi{\usefont{OT1}{cmss}{m}{sl}}
  \DeclareMathAlphabet\mathsfi            {OT1}{cmss}{m}{sl}
  \DeclareTextFontCommand\textsfb{\usefont{OT1}{cmss}{bx}{n}}
  \DeclareMathAlphabet\mathsfb            {OT1}{cmss}{bx}{n}
  \DeclareTextFontCommand\textsfbi{\usefont{OT1}{cmss}{m}{sl}}
  \DeclareMathAlphabet\mathsfbi            {OT1}{cmss}{m}{sl}
\fi

\DeclareMathSymbol{\varGamma}{\mathord}{letters}{"00}
\DeclareMathSymbol{\varDelta}{\mathord}{letters}{"01}
\DeclareMathSymbol{\varTheta}{\mathord}{letters}{"02}
\DeclareMathSymbol{\varLambda}{\mathord}{letters}{"03}
\DeclareMathSymbol{\varXi}{\mathord}{letters}{"04}
\DeclareMathSymbol{\varPi}{\mathord}{letters}{"05}
\DeclareMathSymbol{\varSigma}{\mathord}{letters}{"06}
\DeclareMathSymbol{\varUpsilon}{\mathord}{letters}{"07}
\DeclareMathSymbol{\varPhi}{\mathord}{letters}{"08}
\DeclareMathSymbol{\varPsi}{\mathord}{letters}{"09}
\DeclareMathSymbol{\varOmega}{\mathord}{letters}{"0A}

\setlength\oddsidemargin{2pc}
\setlength\evensidemargin{2pc}
\setlength\marginparwidth{2.0cm}
\setlength\marginparsep{10\p@}

\setlength\topmargin{2.5pc}
\setlength\headheight{12\p@}
\setlength\headsep{8\p@}
\setlength\topskip{12\p@}
\setlength\footskip{\z@}

\setlength\textheight{44\baselineskip}
\addtolength\textheight{\topskip}
\setlength\textwidth{32pc}
\setlength\columnsep{10\p@}
\setlength\columnseprule{\z@}

\setlength\footnotesep{\z@}
\setlength{\skip\footins}{9\p@ \@plus 3\p@}

\setlength\floatsep{10\p@ \@plus 2\p@ \@minus 1\p@}
\setlength\textfloatsep{12\p@ \@plus 6\p@}
\setlength\intextsep{\floatsep}
\setlength\dblfloatsep{\floatsep}
\setlength\dbltextfloatsep{\floatsep}
\setlength\@fptop{\z@ \@plus 0fil}
\setlength\@fpsep{10\p@ \@plus 0fil}
\setlength\@fpbot{\z@ \@plus 3fil}
\setlength\@dblfptop{\z@ \@plus 0fil}
\setlength\@dblfpsep{10\p@ \@plus 0fil}
\setlength\@dblfpbot{\z@ \@plus 3fil}
\setlength\marginparpush{5\p@}

\setlength\parskip{\z@ \@plus .3\p@}
\setlength\parindent{1em}
\setlength\partopsep{\z@ \@plus 1\p@}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\@beginparpenalty -\@lowpenalty
\@endparpenalty -\@lowpenalty
\@itempenalty -\@lowpenalty
\clubpenalty\z@
\widowpenalty\@M

\newcommand\partname{Part}
\newcommand\part{\par\addvspace{4ex}\@afterindentfalse \secdef\@part\@spart}

\def\@part[#1]#2{%
  \ifnum \c@secnumdepth >\m@ne
    \refstepcounter{part}%
    \addcontentsline{toc}{part}{\partname\ \thepart: #1}%
  \else
    \addcontentsline{toc}{part}{#1}%
  \fi
  {\parindent \z@ \centering
   \ifnum \c@secnumdepth >\m@ne
     \normalfont\large\rmfamily \MakeUppercase{\partname}\ %
     \ifcase\thepart \or ONE \or TWO \or THREE \or FOUR \or FIVE
       \or SIX \or SEVEN \or EIGHT \or NINE \or TEN \else \fi
   \par \nobreak
   \fi
   \normalfont\LARGE\rmfamily #2 \markboth{}{}\par}%
  \nobreak \vskip 3ex \@afterheading
}

\def\@spart#1{%
  {\parindent \z@ \centering\normalfont\LARGE\rmfamily #1\par}%
  \nobreak \vskip 3ex \@afterheading
}

\newcommand\@raggedright{\rightskip=\@flushglue \let\\=\@centercr
  \parindent\z@}
\newcommand\@raggedcenter{\centering\rightskip=\@flushglue
  \let\\=\@centercr}

\newif\ifno@eqnreset
\newcommand\nosectioneqnreset{\global\no@eqnresettrue
  \renewcommand\theequation{\@arabic\c@equation}%
}

\newcommand\section{%
  \ifno@eqnreset
    \ifappendix \setcounter{equation}\z@\fi
  \else
    \setcounter{equation}\z@
  \fi
  \@startsection {section}{1}{\z@}
    {-19.5\p@ \@plus -6\p@ \@minus -2\p@}
    {6\p@}
    {\normalfont\normalsize\bfseries\@raggedcenter}%
}

\newcommand\subsection{%
  \@startsection{subsection}{2}{\z@}
    {-19.5\p@ \@plus -3\p@ \@minus -3\p@}
    {6\p@}
    {\normalfont\normalsize\bfseries\@raggedcenter}%
}

\newcommand\subsubsection{%
  \@startsection{subsubsection}{3}{\z@}
    {-19.5\p@ \@plus -3\p@ \@minus -3\p@}
    {6\p@}
    {\normalfont\normalsize\itshape\raggedright}%
}

\newcommand\paragraph{%
  \@startsection{subsubsection}{4}{\z@}
    {-19.5\p@ \@plus -3\p@ \@minus -3\p@}
    {0.1\p@}
    {\normalfont\normalsize\bfseries\raggedright}%
}

\newcommand\subparagraph{}

\def\@seccntformat#1{% FROM LATEX.LTX
  \normalfont
  \ifnum\sec@num=\thr@@ \rmfamily \else \bfseries\fi
  \csname the#1\endcsname\enskip
}

\def\@secappcntformat#1{%
  \ifappendix \appendixname~\fi
  \csname the#1\endcsname\relax
  \ifx\apphe@d\@empty \else \enskip\fi
}

\newcount\sec@num

\def\@sect#1#2#3#4#5#6[#7]#8{% FROM LATEX.LTX
  \sec@num=#2%
  \ifnum #2>\c@secnumdepth
    \let\@svsec\@empty
  \else
    \refstepcounter{#1}%
    \ifnum #2>\@ne
      \protected@edef\@svsec{\@seccntformat{#1}\relax}%
    \else
      \def\apphe@d{#8}%
      \protected@edef\@svsec{\@secappcntformat{#1}\relax}%
    \fi
  \fi
  \@tempskipa #5\relax
  \ifdim \@tempskipa>\z@
    \begingroup
    #6{%
       \@hangfrom{\hskip #3\relax\@svsec}%
         \interlinepenalty \@M #8\@@par}%
    \endgroup
    \csname #1mark\endcsname{#7}%
    \ifnum #2=\@ne
      \addcontentsline{toc}{#1}{%
        \ifnum #2>\c@secnumdepth \else
        \ifappendix \appendixname~\csname the#1\endcsname
        \else \protect\numberline{\csname the#1\endcsname.}#7\fi
        \fi}%
    \else
      \addcontentsline{toc}{#1}{%
        \ifnum #2>\c@secnumdepth \else
          \protect\numberline{\csname the#1\endcsname.}#7\fi}%
    \fi
  \else
    \def\@svsechd{%
      #6{\hskip #3\relax
      \@svsec #8}%
      \csname #1mark\endcsname{#7}%
      \addcontentsline{toc}{#1}{%
        \ifnum #2>\c@secnumdepth \else
          \protect\numberline{\csname the#1\endcsname.}%
        \fi
        #7}}%
  \fi
  \@xsect{#5}}

\newif\ifappendix
\newcommand\appendixname{Appendix}

\newcommand\appendix{\par
  \@addtoreset{figure}{section}%
  \@addtoreset{table}{section}%
  \setcounter{section}\z@
  \renewcommand\thesection{\@Alph\c@section}%
  \renewcommand\theequation{\thesection\,\@arabic\c@equation}%
  \renewcommand\thefigure  {\thesection\,\@arabic\c@figure}%
  \renewcommand\thetable   {\thesection\,\@arabic\c@table}%
  \global\appendixtrue
}

\setcounter{secnumdepth}{3}

\newcounter{part}
\newcounter{section}
\newcounter{subsection}[section]
\newcounter{subsubsection}[subsection]
\renewcommand\thepart          {\@arabic\c@part}
\renewcommand\thesection       {\@arabic\c@section}
\renewcommand\thesubsection    {\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection {\thesubsection.\@arabic\c@subsubsection}

\setlength\leftmargini{\z@}
\setlength\labelsep{\z@}
\setlength\labelwidth{\z@}

\newcommand\makeitemlabel[1]{\hspace\labelwidth\hspace\labelsep\hfill #1}
\newcommand\makeRRlabel[1]{\hss\llap{\normalfont\rmfamily #1}}

\def\@listI{\parsep\z@ \topsep\z@ \partopsep\z@ \itemsep\z@ \rightmargin\z@
 \leftmargin\labelwidth \advance\leftmargin\labelsep}
\def\@listi{\@listI \leftmargin\leftmargini \itemindent 1em}
\def\@listii{\@listI \advance\leftmargin 1em}
\def\@listiii{\@listI}
\def\@listiv{\@listI}
\def\@listv{\@listI}
\def\@listvi{\@listI}
\@listi

\def\itemize{% FROM LATEX.LTX
  \ifnum \@itemdepth >\thr@@ \@toodeep\else
    \advance\@itemdepth \@ne
    \edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
    \expandafter
    \list
      \csname\@itemitem\endcsname
      {\labelsep .5em%
       \settowidth\leftmargin{\csname\@itemitem\endcsname}%
       \itemindent -\leftmargin
       \advance\leftmargin \labelsep
       \listparindent\parindent \topsep 0.5\baselineskip
       \parsep \z@ \@plus .5\p@}%
  \fi
}

\newcommand\labelitemi{$\m@th\bullet$}
\newcommand\labelitemii{$\m@th\circ$}
\newcommand\labelitemiii{\normalfont\bfseries ---}
\newcommand\labelitemiv{\normalfont\bfseries --}

\newcommand\makeRRenumlabel[1]{\hss\llap{{\normalfont\rmfamily #1}}}

\def\enumerate{% FROM LATEX.LTX
  \ifnum \@enumdepth >\thr@@ \@toodeep\else
    \advance\@enumdepth \@ne
    \edef\@enumctr{enum\romannumeral\the\@enumdepth}%
  \fi
  \@ifnextchar [{\@enumeratetwo}{\@enumerateone}%
}

\def\@enumerateone{%
  \expandafter
  \list
    \csname label\@enumctr\endcsname
    {\usecounter{\@enumctr}\leftmargin 2.7em\labelsep 0.5em\itemindent \z@
     \listparindent\parindent \topsep 0.5\baselineskip \parsep 3\p@
     \let\makelabel\makeRRenumlabel}%
}

\def\@enumeratetwo[#1]{%
  \expandafter
  \list
    \csname label\@enumctr\endcsname
    {\settowidth\labelwidth{\normalfont\rmfamily (#1)}%
     \labelsep 0.5em\itemindent \z@
     \leftmargin\labelwidth \advance\leftmargin\labelsep
     \listparindent\parindent \topsep 0.5\baselineskip \parsep 3\p@
     \usecounter{\@enumctr}%
    \let\makelabel\makeRRenumlabel}%
}

\newcommand\labelenumi  {{\normalfont\rmfamily(\theenumi)}}
\newcommand\labelenumii {{\normalfont\rmfamily(\theenumii)}}
\newcommand\labelenumiii{{\normalfont\rmfamily(\theenumiii)}}
\newcommand\labelenumiv {{\normalfont\rmfamily(\theenumiv)}}

\renewcommand\theenumi{\@arabic\c@enumi}
\renewcommand\theenumii{\@alph\c@enumii}
\renewcommand\theenumiii{\@roman\c@enumiii}
\renewcommand\theenumiv{\@Roman\c@enumiv}

\renewcommand\p@enumii{\theenumi.}
\renewcommand\p@enumiii{\p@enumii\theenumii.}
\renewcommand\p@enumiv{\p@enumiii\theenumiii.}

\newcommand*\descriptionlabel[1]{\hspace\labelsep \normalfont\bfseries #1}

\newenvironment{description}
  {\list{}{\leftmargin 2.7em\labelsep 0.5em\itemindent -2.7em%
   \listparindent\parindent \topsep 0.5\baselineskip \parsep 3\p@
  \let\makelabel\descriptionlabel}}
  {\endlist}

\newenvironment{verse}
  {\let\\=\@centercr
   \list{}{\itemsep\z@
   \itemindent -2.5em%
   \listparindent \itemindent
   \rightmargin\leftmargin
   \advance\leftmargin 2.5em}\item[]}
  {\endlist}

\newenvironment{quotation}
  {\list{}{\listparindent\parindent
   \itemindent\z@ \leftmargin\z@ \rightmargin\leftmargin
   \partopsep\z@
   \topsep 13\p@ \@plus 1\p@ \@minus 2\p@
   \parsep \z@}\item[]\strut\quotesize\ignorespaces}
  {\endlist}

\let\quote\quotation
\let\endquote\endquotation

\newif\if@intheoremlike

\def\@xnthm#1#2[#3]{% FROM LATEX.LTX
\expandafter\@ifdefinable\csname #1\endcsname
  {\@definecounter{#1}\@newctr{#1}[#3]%
   \expandafter\xdef\csname the#1\endcsname{\expandafter\noexpand
   \csname the#3\endcsname \@thmcountersep \@thmcounter{#1}}%
   \global\@namedef{#1}{\@thm{#1}{#2}}%
   \global\@namedef{#1*}{\@thmstar{#1}{#2}}%
   \global\@namedef{end#1}{\@endtheorem}%
   \global\@namedef{end#1*}{\@endthmstar}}%
}

\def\@ynthm#1#2{% FROM LATEX.LTX
\expandafter\@ifdefinable\csname #1\endcsname
  {\@definecounter{#1}%
   \expandafter\xdef\csname the#1\endcsname{\@thmcounter{#1}}%
   \global\@namedef{#1}{\@thm{#1}{#2}}%
   \global\@namedef{#1*}{\@thmstar{#1}{#2}}%
   \global\@namedef{end#1}{\@endtheorem}%
   \global\@namedef{end#1*}{\@endthmstar}}%
}

\def\@othm#1[#2]#3{% FROM LATEX.LTX
  \@ifundefined{c@#2}{\@nocounterr{#2}}%
    {\expandafter\@ifdefinable\csname #1\endcsname
    {\global\@namedef{the#1}{\@nameuse{the#2}}%
     \global\@namedef{#1}{\@thm{#2}{#3}}%
     \global\@namedef{#1*}{\@thmstar{#2}{#3}}%
     \global\@namedef{end#1}{\@endtheorem}%
     \global\@namedef{end#1*}{\@endthmstar}}}%
}

\newif\ifrembrks
\newcommand\removebrackets{\rembrkstrue}

\def\thm@body{\topsep 1\baselineskip \@plus 1\p@
  \labelwidth\z@ \itemindent\z@ \labelsep .5em%
  \listparindent\parindent \let\makelabel\makeitemlabel
}

\def\@begintheorem#1#2{% FROM LATEX.LTX
  \list{}{\thm@body}%
  \item[{\normalfont\bfseries #1\ #2}]\normalfont\itshape
}

\def\@opargbegintheorem#1#2#3{% FROM LATEX.LTX
  \list{}{\thm@body}%
  \item[{\normalfont\bfseries #1\ #2}\ {\normalfont\rmfamily
     \ifrembrks #3\global\rembrksfalse\else (#3)\fi}]%
  \normalfont\itshape
}

\def\@endtheorem{\endlist\par\@doendpe}% FROM LATEX.LTX

\def\@thmstar#1#2{\@ifnextchar[{\@ythmstar{#1}{#2}}{\@xthmstar{#1}{#2}}}
\def\@xthmstar#1#2{\@beginthmstar{#2}{\csname the#1\endcsname}\ignorespaces}
\def\@ythmstar#1#2[#3]{\@opargbeginthmstar{#2}{\csname the#1\endcsname}{#3}%
  \ignorespaces}

\def\@thmstarcounter#1{{}}
\def\@thmstarcountersep{{}}

\def\@beginthmstar#1#2{%
  \list{}{\thm@body}%
  \item[{\normalfont\bfseries #1}]\normalfont\itshape
}

\def\@opargbeginthmstar#1#2#3{%
  \list{}{\thm@body}%
  \item[{\normalfont\bfseries #1}\ {\normalfont\rmfamily
    \ifrembrks #3\global\rembrksfalse\else (#3)\fi}]%
  \normalfont\itshape
}

\def\@endthmstar{\endlist\par\@doendpe}

\newcommand\newdefinition{\newnumbered}
\newcommand\newnumbered[1]{\@ifnextchar[{\@oexmpl{#1}}{\@nexmpl{#1}}}
\def\@nexmpl#1#2{\@ifnextchar[{\@xnexmpl{#1}{#2}}{\@ynexmpl{#1}{#2}}}

\def\@xnexmpl#1#2[#3]{%
\expandafter\@ifdefinable\csname #1\endcsname
  {\@definecounter{#1}\@newctr{#1}[#3]%
   \expandafter\xdef\csname the#1\endcsname{\expandafter\noexpand
   \csname the#3\endcsname \@exmplcountersep \@exmplcounter{#1}}%
   \global\@namedef{#1}{\@exmpl{#1}{#2}}%
   \global\@namedef{#1*}{\@exmplstar{#1}{#2}}%
   \global\@namedef{end#1}{\@endexample}%
   \global\@namedef{end#1*}{\@endexmplstar}}%
}

\def\@ynexmpl#1#2{%
\expandafter\@ifdefinable\csname #1\endcsname
  {\@definecounter{#1}%
   \expandafter\xdef\csname the#1\endcsname{\@exmplcounter{#1}}%
   \global\@namedef{#1}{\@exmpl{#1}{#2}}%
   \global\@namedef{#1*}{\@exmplstar{#1}{#2}}%
   \global\@namedef{end#1}{\@endexample}%
   \global\@namedef{end#1*}{\@endexmplstar}}%
}

\def\@oexmpl#1[#2]#3{%
  \@ifundefined{c@#2}{\@nocounterr{#2}}%
    {\expandafter\@ifdefinable\csname #1\endcsname
    {\global\@namedef{the#1}{\@nameuse{the#2}}%
     \global\@namedef{#1}{\@exmpl{#2}{#3}}%
     \global\@namedef{#1*}{\@exmplstar{#2}{#3}}%
     \global\@namedef{end#1}{\@endexample}
     \global\@namedef{end#1*}{\@endexmplstar}}}%
}

\def\@exmpl#1#2{\refstepcounter{#1}%
  \@ifnextchar[{\@yexmpl{#1}{#2}}{\@xexmpl{#1}{#2}}%
}

\def\@xexmpl#1#2{\@beginexample{#2}{\csname the#1\endcsname}\ignorespaces}
\def\@yexmpl#1#2[#3]{\@opargbeginexample{#2}{\csname the#1\endcsname}{#3}%
  \ignorespaces}

\def\@beginexample#1#2{%
  \list{}{\thm@body}%
  \item[{\normalfont\bfseries #1\ #2}]\normalfont\rmfamily
}

\def\@opargbeginexample#1#2#3{%
  \list{}{\thm@body}%
  \item[{\normalfont\bfseries #1\ #2}\ {\normalfont\rmfamily
    \ifrembrks #3\global\rembrksfalse\else (#3)\fi}]%
  \normalfont\rmfamily
}

\def\@endexample{\endlist\par\@doendpe}
\def\@exmplcounter#1{\noexpand\arabic{#1}}
\def\@exmplcountersep{.}

\def\@exmplstar#1#2{\@ifnextchar[{\@yexmplstar{#1}{#2}}{\@xexmplstar{#1}{#2}}}
\def\@xexmplstar#1#2{\@beginexmplstar{#2}{\csname the#1\endcsname}%
  \ignorespaces}

\def\@yexmplstar#1#2[#3]{\@opargbeginexmplstar{#2}{\csname the#1\endcsname}%
  {#3}\ignorespaces}

\def\@exmplstarcounter#1{{}}
\def\@exmplstarcountersep{{}}

\def\@beginexmplstar#1#2{%
  \list{}{\thm@body}%
  \item[{\normalfont\bfseries #1}]\normalfont\rmfamily
}

\def\@opargbeginexmplstar#1#2#3{%
  \list{}{\thm@body}%
  \item[{\normalfont\bfseries #1}\ {\normalfont\rmfamily
         \ifrembrks #3\global\rembrksfalse\else (#3)\fi}]%
  \normalfont\rmfamily
}

\def\@endexmplstar{\endlist\par\@doendpe}

\newsavebox{\proofsavebox}

\ifprodtf
  \sbox{\proofsavebox}{$\CUPproofbox$}
  \newcommand\proofbox{\hbox{$\CUPproofbox$}}
\else
  \sbox{\proofsavebox}
    {\unitlength 1pt\begin{picture}(6.5,6.5)%
     \put(0,0){\framebox(6.5,6.5){}}\end{picture}}
  \newcommand\proofbox{\usebox{\proofsavebox}\relax}
\fi

\newlength\rpbdim
\newcommand\raiseproofboxby[1]{\setlength\rpbdim{#1}}
\def\raiseproofb@x{\vbox to\z@{\hsize=\wd\proofsavebox\vss
  \hbox{\proofbox}\vskip\rpbdim}%
  \setlength\rpbdim{\z@}%
}

\newenvironment{proof}
  {\@ifnextchar[{\proofwith}{\proofwithout}}
  {\hspace*{1em \@plus 1000fil}\nobreak
   \ifdim\rpbdim=\z@ \hbox{\proofbox}\par \vskip\baselineskip
   \else \raiseproofb@x\par\fi \@doendpe}

\newenvironment{proof*}
  {\proof}
  {\par\vskip\baselineskip\@doendpe}

\def\proofwith[#1]{\par\addvspace{\@ne\baselineskip \@plus 1\p@}%
   \noindent {\normalfont\bfseries #1}\enskip
   \ignorespaces}

\def\proofwithout{\par\addvspace{\@ne\baselineskip \@plus 1\p@}%
   \noindent {\normalfont\bfseries Proof}\enskip
   \ignorespaces}

\renewcommand\theequation{\thesection.\@arabic\c@equation}

% subeqn.sty follows
% Copyright (c) 1990 by Stephen Gildea
% Permission to copy all or part of this work is granted, provided
% that the copies are not made or distributed for resale, and that
% the copyright notice and this notice are retained.
% Modified for use with EJM/LaTeX2e.

\newif\ifuseAMSsubequ
\newcommand\useAMSsubequations{\global\useAMSsubequtrue}

\newtoks\@stequation

\newenvironment{CUPsubequations}
  {\refstepcounter{equation}%
   \edef\@savedequation{\the\c@equation}%
   \@stequation=\expandafter{\theequation}%
   \edef\@savedtheequation{\the\@stequation}%
   \edef\oldtheequation{\theequation}%
   \setcounter{equation}\z@
   \def\theequation{\oldtheequation\,{\itshape\@alph\c@equation}}}
  {\setcounter{equation}{\@savedequation}%
   \@stequation=\expandafter{\@savedtheequation}%
   \edef\theequation{\the\@stequation}\global\@ignoretrue}

\newcommand\returnthesubequation{%
  \gdef\thesubequation{\theequation\,{\itshape\@alph\c@subequation}}}

% subeqnarray.sty
% Copyright (C) 1988 by Johannes Braams, Dr Neher Laboratories
% Modified for use with EJM/LaTeX2e.

\newcounter{subequation}[equation]
\renewcommand\thesubequation{\theequation\,{\itshape\@alph\c@subequation}}
\def\@subeqnnum{{\normalfont\rmfamily (\thesubequation)}}

\newcommand\slabel[1]{\@bsphack % MODIFIED IN LINE WITH LATEX.LTX
  \protected@write\@auxout{}%
    {\string\newlabel{#1}{{\thesubequation}{\thepage}}}%
  \@esphack}

\newenvironment{subeqnarray} % MODIFIED IN LINE WITH LATEX.LTX
  {\stepcounter{equation}%
   \let\@currentlabel=\theequation
   \global\c@subequation\@ne
   \global\@eqnswtrue
   \m@th
   \global\@eqcnt\z@
   \tabskip\@centering
   \let\\\@subeqncr
   $$\everycr{}\halign to\displaywidth\bgroup
     \hskip\@centering$\displaystyle\tabskip\z@skip{##}$\@eqnsel
    &\global\@eqcnt\@ne \hfil${\;##\;}$\hfil
    &\global\@eqcnt\tw@
     $\displaystyle\tabskip\z@skip{##}$\hfil\tabskip\@centering
    &\hb@xt@\z@\bgroup\hss##\egroup
      \tabskip\z@skip \cr}
  {\@@subeqncr\egroup $$\global\@ignoretrue}

\def\@subeqncr{{\ifnum0=`}\fi\@ifstar{\global\@eqpen\@M
    \@ysubeqncr}{\global\@eqpen\interdisplaylinepenalty \@ysubeqncr}}

\def\@ysubeqncr{\@ifnextchar [{\@xsubeqncr}{\@xsubeqncr[\z@]}}

\def\@xsubeqncr[#1]{\ifnum0=`{\fi}\@@subeqncr
   \noalign{\penalty\@eqpen\vskip\jot\vskip #1\relax}}

\def\@@subeqncr{\let\@tempa\relax
    \ifcase\@eqcnt \def\@tempa{& & &}\or \def\@tempa{& &}
      \else \def\@tempa{&}\fi
     \@tempa \if@eqnsw\@subeqnnum\refstepcounter{subequation}\fi
     \global\@eqnswtrue\global\@eqcnt\z@\cr}

\let\@ssubeqncr=\@subeqncr
\newenvironment{subeqnarray*}
  {\def\@subeqncr{\nonumber\@ssubeqncr}\subeqnarray}
  {\global\advance\c@equation\m@ne \nonumber\endsubeqnarray}

\renewenvironment{eqnarray} % MODIFIED IN LINE WITH LATEX.LTX
  {\stepcounter{equation}%
   \let\@currentlabel=\theequation
   \global\@eqnswtrue
   \m@th
   \global\@eqcnt\z@
   \tabskip\@centering
   \let\\\@eqncr
   $$\everycr{}\halign to \displaywidth\bgroup
     \hskip\@centering$\displaystyle\tabskip\z@skip{##}$\@eqnsel
   &\global\@eqcnt\@ne \hfil${\;##\;}$\hfil
   &\global\@eqcnt\tw@
     $\displaystyle\tabskip\z@skip{##}$\hfil \tabskip\@centering
   &\hb@xt@\z@\bgroup\hss##\egroup
     \tabskip\z@skip \cr}
  {\@@eqncr\egroup\global\advance\c@equation\m@ne
   $$\global\@ignoretrue}

\setlength\arraycolsep{5\p@}
\setlength\tabcolsep{3\p@}
\setlength\arrayrulewidth{.25\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\setlength{\skip\@mpfootins}{5\p@ \@plus .4\p@ \@minus .4\p@}
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.5\p@}

\newcommand\affiliation[1]{\gdef\@affiliation{#1}}
\gdef\@affiliation{}

\newcommand\maketitle{\par
  \begingroup
    \newpage
    \global\@topnum\z@
    \def\thefootnote{\mbox{$\@fnsymbol\c@footnote$}}%
    \let\@makefntext\tp@makefntext
    \def\@makefnmark{\hbox{$\@thefnmark$}}%
    \@maketitle
    \thispagestyle{titlepage}\@thanks
  \endgroup
  \setcounter{footnote}\z@
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@title\@empty
  \global\let\@author\@empty
  \global\let\@affiliation\@empty
  \global\let\maketitle\relax
  \global\let\thanks\relax
}

\renewcommand\and{{\normalfont\normalsize\rmfamily and\ns}}

\def\@maketitle{%
  \newpage
  \vspace*{10\p@}%
  {\centering \sloppy
   {\normalfont\LARGE\bfseries \strut\@title \par}%
   \vskip 13\p@ \@plus 2\p@ \@minus 1\p@
   {\normalfont\normalsize\rmfamily \strut\@author \par}%
   \vskip 10\p@ \@plus 1\p@
   {\normalfont\authorsize\itshape \strut\@affiliation \par}%
   \vskip 17\p@ \@plus 2\p@ \@minus 1\p@
   {\normalfont\authorsize\rmfamily \strut
     ({\normalfont\itshape Received \@date})}%
  \par}%
  \vskip 19\p@ \@plus 2\p@ \@minus 1\p@
}

\newenvironment{abstract}
  {\par\leftskip 1pc\rightskip\leftskip
   \normalfont\abstractsize\noindent\ignorespaces}
  {\par\vskip 9\p@ \@plus 1\p@ \@minus 1\p@}

\newenvironment{keywords}
  {\par\noindent{\bf Key Words:}\ \normalfont\abstractsize\noindent\ignorespaces}
  {\par\vskip 9\p@ \@plus 1\p@ \@minus 1\p@
   \vbox{\centerline{\rule[4\p@]{30pc}{.4\p@}}}}

\mark{{}{}}
\renewcommand\author{\@ifnextchar [{\@authortwo}{\@authorone}}
\def\@authortwo[#1]#2{\gdef\@author{#2}\gdef\@shortauthor{#1}}
\def\@authorone#1{\gdef\@author{#1}%
  \begingroup\def\and{and}\gdef\@shortauthor{#1}\endgroup}
\gdef\@author{\mbox{}}

\newcommand\shortauthor[1]{\gdef\@shortauthor{#1}}
\gdef\@shortauthor{}

\renewcommand\title{\@ifnextchar [{\@titletwo}{\@titleone}}
\def\@titletwo[#1]#2{\gdef\@title{#2}\gdef\@shorttitle{#1}}
\def\@titleone#1{\gdef\@title{#1}\gdef\@shorttitle{#1}}
\gdef\@title{\mbox{}}

\newcommand\shorttitle[1]{\gdef\@shorttitle{#1}}
\gdef\@shorttitle{}

\newcommand\volume[1]{\gdef\@volume{#1}}
\gdef\@volume{00}

\newcommand\pagerange[1]{\gdef\@pagerange{#1}}
\gdef\@pagerange{1--000}

\newcommand\pubyear[1]{\def\@pubyear{#1}}
\gdef\@pubyear{19XX}

\newcommand\doi[1]{\gdef\@doi{10.1017/#1}}
\gdef\@doi{10.1017/S000000000000000}

\newcommand\journal[1]{\gdef\@journal{#1}}
\def\@journal{%
  \vbox to 5.5\p@{\noindent\parbox[t]{5in}{\normalfont\catchlinesize
    {\itshape Euro.\ Jnl of Applied Mathematics} (\@pubyear), {\itshape vol.}\ %
    \@volume, {\itshape pp.} \@pagerange.\quad \copyright\ \@pubyear\ 
    Cambridge University Press\\[2\p@]
    DOI: \@doi\quad {Printed in the United Kingdom}}%
  \vss}%
}

\gdef\@underjournal{%
  \vbox to 5.5\p@{\noindent\parbox[t]{5in}{\normalfont\catchlinesize
    {\itshape Under consideration for publication in
     Euro.\ Jnl of Applied Mathematics}\\[2.5\p@]
      {\ \ }}%
  \vss}%
}

\def\ps@headings{\let\@mkboth\markboth
  \def\@oddhead{\hfil{\itshape\@shorttitle}\hfil\llap{\thepage}}%
  \def\@evenhead{\rlap{\thepage}\hfil\itshape\@shortauthor\hfil}%
  \let\@oddfoot\@empty
  \let\@evenfoot\@oddfoot
  \def\sectionmark##1{\markboth{##1}{}}%
  \def\subsectionmark##1{\markright{##1}}%
}

\def\ps@myheadings{\let\@mkboth\@gobbletwo
  \def\@oddhead{\hfil{\itshape\rightmark}\hfil\llap{\thepage}}%
  \def\@evenhead{\rlap{\thepage}\hfil\itshape\leftmark\hfil}%
  \let\@oddfoot\@empty
  \let\@evenfoot\@oddfoot
  \let\sectionmark\@gobble
  \let\subsectionmark\@gobble
}

\ifprodtf
  \def\@j@urnal{\@journal}
\else
  \def\@j@urnal{\@underjournal}
\fi

\def\ps@titlepage{\leftskip\z@\let\@mkboth\@gobbletwo
  \def\@oddhead{\@j@urnal \hfil\llap{\thepage}}%
  \let\@evenhead\@oddhead
  \let\@oddfoot\@empty
  \let\@evenfoot\@oddfoot
  \let\sectionmark\@gobble
  \let\subsectionmark\@gobble
}

\newcommand\@pnumwidth{1.55em}
\newcommand\@tocrmarg {2.55em}
\newcommand\@dotsep{4.5}
\setcounter{tocdepth}{2}

\def\@dottedtocline#1#2#3#4#5{% FROM LATEX.LTX
  \ifnum #1>\c@tocdepth \else
    \vskip \z@ \@plus.2\p@
    {\leftskip #2\relax \rightskip \@tocrmarg \parfillskip -\rightskip
     \parindent #2\relax\@afterindenttrue
     \interlinepenalty\@M
     \leavevmode
     \@tempdima #3\relax
     \advance\leftskip \@tempdima \null\hskip -\leftskip
     {#4}\nobreak
     \leaders\hbox{$\m@th
        \mkern \@dotsep mu\hbox{\phantom{.}}\mkern \@dotsep
        mu$}\hfill
     \nobreak
     \hb@xt@\@pnumwidth{\hfil\normalfont \normalcolor #5}%
     \par}%
  \fi}

\newcommand\contentsname{Contents}

\newcommand\tableofcontents{%
  \@restonecolfalse
  \if@twocolumn\@restonecoltrue\onecolumn\fi
  \section*{\contentsname}%
  \@starttoc{toc}%
  \if@restonecol\twocolumn\fi \par\vspace{12\p@}%
}

\newcommand*\l@section[2]{%
  \ifnum \c@tocdepth >\z@
    \addpenalty{\@secpenalty}%
    \@tempdima 1.5em%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth \normalfont\rmfamily \leavevmode
      \advance\leftskip\@tempdima \hskip -\leftskip
      {\normalfont\bfseries #1}\nobreak\hfil
      \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
    \endgroup
  \fi
}

\newcommand*\l@subsection{\@dottedtocline{2}{2.5em}{2.3em}}
\newcommand*\l@subsubsection{\@dottedtocline{3}{3.8em}{3.2em}}

\newif\if@restonecol
\newcommand\indexname{Index}

\newenvironment{theindex}
  {\section*{\indexname}%
   \addcontentsline{toc}{section}{\indexname}%
   \normalfont\footnotesize \parindent\z@ \parskip\z@ \@plus .1\p@\relax
   \let\item\@idxitem}
  {\if@restonecol\onecolumn\else\clearpage\fi}

\newcommand\@idxitem{\par\hangindent 1em\relax}
\newcommand\subitem{\@idxitem\hspace*{1em}}
\newcommand\subsubitem{\subitem}
\newcommand\indexspace{\par\vskip\baselineskip \relax}

\renewcommand\footnoterule{\kern-3\p@ \hrule \@width\z@ \kern 2.6\p@}

\long\def\@mpfootnotetext#1{% FROM LATEX.LTX
  \global\setbox\@mpfootins\vbox{%
    \unvbox\@mpfootins
    \reset@font\footnotesize
    \hsize\columnwidth
    \@parboxrestore
    \protected@edef\@currentlabel
         {\csname p@mpfootnote\endcsname\@thefnmark}%
    \color@begingroup
      \@makempfntext{%
        \rule\z@\footnotesep\ignorespaces#1\@finalstrut\strutbox}%
    \color@endgroup}}

\newcommand\tp@makefntext[1]{%
  \strut\hb@xt@ 1.5em{\hss$\@thefnmark$}\enskip \normalfont\rmfamily #1}
\newcommand\@makefntext[1]{%
  \strut\hb@xt@ 1.5em{\hss$^{\@thefnmark}$}\enskip \normalfont\rmfamily #1}
\newcommand\@makempfntext[1]{%
  \setbox\@tempboxa\hbox{#1}%
  \ifdim\wd\@tempboxa>\textwidth
    \strut\hb@xt@ 1.5em{\hss$^{\@thefnmark}$}%
  \else
    \centering $^{\@thefnmark}$%
  \fi
  \enskip \normalfont\rmfamily #1%
}

\renewcommand\@makefnmark{\hbox{$^{\@thefnmark}$}}

\renewcommand\thefootnote{\@arabic\c@footnote}% FROM LATEX.LTX
\renewcommand\thempfootnote{\mathrm{\@alph\c@mpfootnote}}% FROM LATEX.LTX

% FROM LATEX.LTX
\def\@fnsymbol#1{\ensuremath{\ifcase#1\or \dagger\or \ddagger\or
  \mathchar "278\or \|\or \mathchar "27B\or \dagger\dagger\or
  \ddagger\ddagger\or \mathchar "278\mathchar "278\or \|\|\or
  \mathchar "27B\mathchar "27B\or \else\@ctrerr\fi}%
}

\setcounter{topnumber}{3}
\renewcommand\topfraction{.9}
\setcounter{bottomnumber}{0}
\renewcommand\bottomfraction{.9}
\setcounter{totalnumber}{3}
\renewcommand\textfraction{.1}
\renewcommand\floatpagefraction{.9}
\setcounter{dbltopnumber}{3}
\renewcommand\dbltopfraction{.9}
\renewcommand\dblfloatpagefraction{.5}

\newif\ifrem@fullpt
\newcommand\removefullpoint{\global\rem@fullpttrue}

\newif\ifbot@fig
\newenvironment{bottomfigure}{\def\fps@figure{b}%
  \setcounter{bottomnumber}{1}%
  \global\bot@figtrue
  \let\@makecaption=\@make@figure@caption
  \@float{figure}\fstyle@figure}
  {\end@float}

\newbox\@tempboxb

\long\def\@make@table@caption#1#2{%
  \begingroup
  \setbox\@tempboxb\hbox{#2}%
  \def\@xtra{\ifdim\wd\@tempboxb>\z@ \ifrem@fullpt\else .\enskip\fi\fi}%
  \setbox\@tempboxa\hbox{#1\@xtra #2}%
  \ifdim\wd\@tempboxa>\tw@\textwidth
    {\let\centering\relax #1\@xtra #2\par}%
  \else
    #1\@xtra #2\par
  \fi
  \endgroup
  \global\rem@fullptfalse
  \vspace{-3.5\p@}%
}

\newcounter{table}
\renewcommand\thetable{\@arabic\c@table}
\def\fps@table{tp}
\def\ftype@table{1}
\def\ext@table{lot}
\newcommand\tablename{Table}
\def\fnum@table{\tablename~\thetable}

\newenvironment{table}
  {\let\@makecaption=\@make@table@caption
   \@float{table}\fstyle@table}
  {\end@float}

\newenvironment{table*}
  {\let\@makecaption=\@make@table@caption
   \@dblfloat{table}\fstyle@table}
  {\end@dblfloat}

\def\fstyle@table{\normalfont\tablebdysize\rmfamily}
\def\fjust@table{\centering}
\def\fcapjust@table{\centering}
\def\fcapsize@table{\normalfont\normalsize\rmfamily}
\def\fcapstyle@table{\normalfont\normalsize\itshape}

\long\def\@make@figure@caption#1#2{%
  \ifbot@fig \rule{\textwidth}{.25\p@}\fi
  \vskip 5\p@ \@plus .4\p@ \@minus .4\p@
  \begingroup
  \setbox\@tempboxb\hbox{#2}%
  \def\@xtra{\ifdim\wd\@tempboxb>\z@ \ifrem@fullpt\else .\enskip\fi\fi}%
  \setbox\@tempboxa\hbox{#1\@xtra #2}%
  \ifdim\wd\@tempboxa>\tw@\textwidth
    {\let\centering\relax #1\@xtra #2\par}%
  \else
    #1\@xtra #2\par
  \fi
  \endgroup
  \global\bot@figfalse
  \global\rem@fullptfalse
}

\newcounter{figure}
\renewcommand\thefigure{\@arabic\c@figure}
\def\fps@figure{tp}
\def\ftype@figure{2}
\def\ext@figure{lof}
\newcommand\figurename{Figure}
\def\fnum@figure{\figurename~\thefigure}

\newenvironment{figure}
  {\let\@makecaption=\@make@figure@caption
   \@float{figure}}
  {\end@float}

\newenvironment{figure*}
  {\let\@makecaption=\@make@figure@caption
   \@dblfloat{figure}}
  {\end@dblfloat}

\def\fstyle@figure{\normalfont\normalsize\rmfamily}
\def\fjust@figure{\centering}
\def\fcapjust@figure{\centering}
\def\fcapsize@figure{\normalfont\figcapsize\scshape}
\def\fcapstyle@figure{\normalfont\figcapsize\rmfamily}

\long\def\@caption#1[#2]#3{% FROM LATEX.LTX
  \par
  \addcontentsline{\csname ext@#1\endcsname}{#1}%
    {\protect\numberline{\csname the#1\endcsname}{\ignorespaces #2}}%
  \begingroup
    \@parboxrestore
    \if@minipage
      \@setminipage
    \fi
    \normalsize
    \@makecaption{\csname fcapjust@#1\endcsname
                  \csname fcapsize@#1\endcsname
                  \csname fnum@#1\endcsname}%
                 {\csname fcapstyle@#1\endcsname \ignorespaces #3}\par
  \endgroup
}

\def\@xfloat #1[#2]{% FROM LATEX.LTX
  \@nodocument
  \def\@captype {#1}%
  \def\@fps {#2}%
  \@onelevel@sanitize \@fps
  \def \reserved@b {!}%
  \ifx \reserved@b \@fps
    \@fpsadddefault
  \else
    \ifx \@fps \@empty
      \@fpsadddefault
    \fi
  \fi
  \ifhmode
    \@bsphack
    \@floatpenalty-\@Mii
  \else
    \@floatpenalty-\@Miii
  \fi
  \ifinner
    \@parmoderr\@floatpenalty\z@
  \else
    \@next\@currbox\@freelist
      {%
       \@tempcnta \sixt@@n
       \expandafter \@tfor \expandafter \reserved@a
         \expandafter :\expandafter =\@fps
         \do
          {%
           \if \reserved@a h%
             \ifodd \@tempcnta
             \else
               \advance \@tempcnta \@ne
             \fi
           \fi
           \if \reserved@a t%
             \@setfpsbit \tw@
           \fi
           \if \reserved@a b%
             \@setfpsbit 4%
           \fi
           \if \reserved@a p%
             \@setfpsbit 8%
           \fi
           \if \reserved@a !%
             \ifnum \@tempcnta>15
               \advance\@tempcnta -\sixt@@n\relax
             \fi
           \fi
           }%
       \@tempcntb \csname ftype@\@captype \endcsname
       \multiply \@tempcntb \@xxxii
       \advance \@tempcnta \@tempcntb
       \global \count\@currbox \@tempcnta
       }%
    \@fltovf
  \fi
  \global \setbox\@currbox
    \color@vbox
      \normalcolor
      \vbox \bgroup
        \hsize\columnwidth
        \@parboxrestore
        \@floatboxreset
        \csname fstyle@\@captype\endcsname
        \csname fjust@\@captype\endcsname
}

\let\oldtabular\tabular
\let\endoldtabular\endtabular

\def\tabular{% FROM LATEX.LTX
  \def\@halignto{to \textwidth}\tabskip\tabcolsep \@plus 1fil\@ttabular
}

\def\@ttabular{\leavevmode \hbox \bgroup $\let\@acol\@tabacol
   \let\@classz\@tabclassz
   \let\@classiv\@tabclassiv \let\\\@tabularcr\@ttabarray}

\def\@ttabarray{\m@th\@ifnextchar[\@tarray{\@ttarray}}

\def\@tarray[#1]#2{\t@barray[#1]{@{\tabskip\tw@\tabcolsep \@plus 3\p@}#2}}
\def\@ttarray#1{\t@barray[c]{@{\tabskip\tw@\tabcolsep \@plus 3\p@}#1}}

\def\t@barray[#1]#2{%
  \if #1t\vtop \else \if#1b\vbox \else \vcenter \fi\fi
  \bgroup
  \setbox\@arstrutbox\hbox{%
    \vrule \@height\arraystretch\ht\strutbox
           \@depth\arraystretch \dp\strutbox
           \@width\z@}%
  \@mkpream{#2}%
  \edef\@preamble{%
    \halign \noexpand\@halignto
      \bgroup \tabskip\z@skip \@arstrut \@preamble \tabskip\tabcolsep \@plus 1fil\cr}%
  \let\@startpbox\@@startpbox \let\@endpbox\@@endpbox
  \let\tabularnewline\\%
  \let\par\@empty
  \let\@sharp##%
  \set@typeset@protect
  \lineskip\z@skip\baselineskip\z@skip
  \@preamble}

\newcommand\ls{\kern.1em\relax}
\newcommand\ns{\kern.4em\relax}

\def\hline{% FROM LATEX.LTX
  \noalign{\ifnum0=`}\fi \vskip 6\p@
  \hrule \@height \arrayrulewidth \vskip 6\p@
  \futurelet \reserved@a\@xhline}

\def\@xhline{% FROM LATEX.LTX
  \ifx\reserved@a\hline
    \vskip -12\p@
    \vskip\doublerulesep
  \fi
  \ifnum0=`{\fi}}

\newenvironment{acknowledgement}{\section*{Acknowledgement}}{\par}
\newenvironment{acknowledgements}{\section*{Acknowledgements}}{\par}

\let\acknowledgments\acknowledgements
\let\endacknowledgments\endacknowledgements
\let\acknowledgment\acknowledgement
\let\endacknowledgment\endacknowledgement

\newcommand\today{}
\edef\today{\number\day\ \ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December
 \fi \ \number\year}

\newif\ifnobr@cket
\DeclareRobustCommand\nbcite{\nobr@ckettrue\cite}

\def\@citex[#1]#2{% FROM LATEX.LTX
  \let\@citea\@empty
  \@cite{\@for\@citeb:=#2\do
    {\@citea\def\@citea{{\normalfont\rmfamily ,}\penalty\@m\ }%
     \edef\@citeb{\expandafter\@firstofone\@citeb}%
     \if@filesw\immediate\write\@auxout{\string\citation{\@citeb}}\fi
     \@ifundefined{b@\@citeb}{\mbox{\reset@font\bfseries ?}%
       \G@refundefinedtrue
       \@latex@warning
         {Citation `\@citeb' on page \thepage \space undefined}}%
       {{\csname b@\@citeb\endcsname}}}}{#1}}

\def\@cite#1#2{% FROM LATEX.LTX
  \ifnobr@cket\else [\fi
  {\normalfont\rmfamily #1\if@tempswa \normalfont\rmfamily, #2\fi}%
  \ifnobr@cket\nobr@cketfalse\else ]\fi
}

\newcommand\refname{References}

\def\@biblabel#1{[{\normalfont\rmfamily #1}]}% FROM LATEX.LTX

\newenvironment{thebibliography}[1]
  {\section*{\refname}%
   \normalfont\bibliosize\rmfamily
   \addcontentsline{toc}{section}{\refname}%
   \list{\@biblabel{\@arabic\c@enumiv}}
     {\settowidth\labelwidth{\@biblabel{#1}}%
      \labelsep 0.5em%
      \leftmargin \labelwidth \advance\leftmargin \labelsep
      \advance\leftmargin 1em%
      \itemsep \z@
      \itemindent 0.5em%
      \advance\itemindent -\leftmargin
      \advance\itemindent \labelwidth
      \usecounter{enumiv}%
      \let\p@enumiv\@empty
      \renewcommand\theenumiv{\@arabic\c@enumiv}%
      \let\makelabel\makeRRlabel}%
   \parindent\z@
   \parskip\z@ \@plus .1\p@
   \sloppy\clubpenalty\z@ \widowpenalty\@M
   \sfcode`\.\@m\relax}
  {\def\@noitemerr
   {\@latex@warning{Empty `thebibliography' environment}}%
   \endlist}

\newcommand\newblock{\hskip .11em \@plus .33em \@minus .07em\relax}

\newif\iffontfound
\newcommand\checkfont[1]{%
  \batchmode
  \font\test=#1\relax
  \errorstopmode
  \fontfoundfalse
  \ifx\test\nullfont \else \fontfoundtrue\fi
}

\pagestyle{headings}
\pagenumbering{arabic}
\frenchspacing
\ifprodtf
  \setlength\oddsidemargin\z@
  \setlength\evensidemargin\z@
  \setlength\topmargin\z@
\fi
\ifreferee
  \renewcommand\baselinestretch{2}\raggedbottom
\else
  \flushbottom
\fi

\endinput

% end of file ejm.cls
