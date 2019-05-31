Kurzanleitung zur Nutzung der Enthaltenen Dateien 
zu Schwerhörendensimulation

Autor: R. Bennett, 2011

---------------------------------------------------------------------------
Anleitung zum Umgang mit den Dateien in diesem Ordner
---------------------------------------------------------------------------
Dieser Ordner enthält Skripte und Funktionen, zum Ausführen einer 
Schwerhörenden- sowie Normalhörendensimulation, wie diese auch
zur Masterarbeit (R. Bennett 'Simulation und Auralisation verminderter 
Frequenzselektivität bei Innenohrschwerhörigkeit') durchgeführt wurden. 

Die enthaltenen Funktionen der Gammatonfilterbank wurden im Rahmen der o.g.
Arbeit von V. Hohmann zur Verfügung gestellt.

---------------------------------------------------------------------------
Anmerkung zu Playrec (bitte beachten!)
---------------------------------------------------------------------------
Die Soundausgabe wurde im vorliegenden Code (Demos) mit der MEX-Datei 
playrec (www.playrec.co.uk) vollzogen. Sollte Playrec nicht auf dem 
Zielsystem installiert sein, oder installiert werden sollen so können 
die entsprechenden Zeilen im Quelltext auskommentiert und ggf. mit 'soundsc' 
ersetzt werden. 

---------------------------------------------------------------------------
Start der Simulation
---------------------------------------------------------------------------
1) mit 'run_Iterations_over1Signal'
Wird die Simulation mit einem Signal gestartet. Am Ende wird ein Plot 
ausgegeben auf dem der Verlauf der Korrelationskoeffizienten über die in 
'f_parames.m' gewählten Iterationen ausgegeben wird. Dabei können alle 
Iterationen angezeigt werden, wenn 'parameters.showAllIterations = 1' 
gesetzt wird. Nur die Baseline, erste und letzte Iteration angezeigt wenn
'parameters.showAllIterations = 0' ist. Es werden zudem die 1. Iteration 
und letzte Iteration der verarbeiten Signale als .wav-File abgespeichert, so 
dass sie auralisiert werden können. (Simulationsergebnisse) 

2) mit 'run_Iterations_overXSignals'
wird der arithmetische Mittelwert der Korralationskoeffizienten über 
X Signale, sowie deren Standardabweichung berechnet. Es werden vier 
Diagramme ausgegeben. In Figure 1 sind Baseline, erste Iteration und 
zehnte Iteration geplottet. Die anderen Diagramme stellen jeweils einen 
Korrelationsverlauf mit den dazugehörigen Standardabweichungen pro Frequenz-
kanal dar. Die verwendeten Signale müssen die gleichen Dateinamen aufweisen
und dieser fortlaufend nummeriert sein (z.B. Sprache1, Sprache2, Sprache3, 
etc.). Die gewünschte Dateibezeichnung kann im Skript
'run_Iterations_overXSignals' verändert werden. 

---------------------------------------------------------------------------
Modell und Algorithmus
---------------------------------------------------------------------------
'f_process.m' ist der Kern der Implementierung und setzt den Algorithmus des 
Modells um. Hier werden die Filterbankobjekte der NH- und 
HI-Analysefilterbanken sowie der NH- und HI-Modelle erzeugt und es findet 
die iterative Berechnung des Korrelationskoeffizienten zwischen den 
Einhüllenden statt.    

---------------------------------------------------------------------------
PATH Umgebung
---------------------------------------------------------------------------
Der Ordner 'CodeDemo' und alle seine Unterordner sind im Matlab Pfad 
einzubinden.

---------------------------------------------------------------------------
Parameter
---------------------------------------------------------------------------
Simulationsparameter werden in der Datei 'f_params.m' den entsprechenden 
Simulationsanforderungen angepasst. Näheres zu den Parametern ist 
in 'f_params' kommentiert.

---------------------------------------------------------------------------
Signals (Sinustöne, Tonkomplexe)
---------------------------------------------------------------------------
Sollte mit Sinustönen oder Tonkomplexen gearbeitet werden, können diese in 
den entsprechenden Zeilen der Datei 'Signals.m' konfiguriert werden.
Alle Einträge in 'Signals.m', die mit Sprache, Störgeräusch und Musik 
zusammenhängen sollen NICHT verändert werden.

---------------------------------------------------------------------------
Audio Daten
---------------------------------------------------------------------------
Der Unterordner AudioData enthält die Audiodateien für die Simulation.
Soundfiles (einkanalig bzw. mono) können hier abgelegt werden.

---------------------------------------------------------------------------
Formatierung von Diagrammen
---------------------------------------------------------------------------
Der Unterordner 'GraphicsFormatting' enthält kleine Funktionen mit denen 
die Diagramme am Ende einer Simulation aufbereitet werden. Diese Funktionen 
werden von den Simulationsskripten aufgerufen.

---------------------------------------------------------------------------
Sound_IO
---------------------------------------------------------------------------
Der Unterordner Sound_IO enthält kleine Funktionen mit denen die Diagramme 
am Ende einer Simulation aufbereitet werden. Diese Funktionen werden von 
den Simulationsskripten aufgerufen.
