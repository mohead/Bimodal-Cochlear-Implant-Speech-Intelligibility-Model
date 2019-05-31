Kurzanleitung zur Nutzung der Enthaltenen Dateien 
zu Schwerh�rendensimulation

Autor: R. Bennett, 2011

---------------------------------------------------------------------------
Anleitung zum Umgang mit den Dateien in diesem Ordner
---------------------------------------------------------------------------
Dieser Ordner enth�lt Skripte und Funktionen, zum Ausf�hren einer 
Schwerh�renden- sowie Normalh�rendensimulation, wie diese auch
zur Masterarbeit (R. Bennett 'Simulation und Auralisation verminderter 
Frequenzselektivit�t bei Innenohrschwerh�rigkeit') durchgef�hrt wurden. 

Die enthaltenen Funktionen der Gammatonfilterbank wurden im Rahmen der o.g.
Arbeit von V. Hohmann zur Verf�gung gestellt.

---------------------------------------------------------------------------
Anmerkung zu Playrec (bitte beachten!)
---------------------------------------------------------------------------
Die Soundausgabe wurde im vorliegenden Code (Demos) mit der MEX-Datei 
playrec (www.playrec.co.uk) vollzogen. Sollte Playrec nicht auf dem 
Zielsystem installiert sein, oder installiert werden sollen so k�nnen 
die entsprechenden Zeilen im Quelltext auskommentiert und ggf. mit 'soundsc' 
ersetzt werden. 

---------------------------------------------------------------------------
Start der Simulation
---------------------------------------------------------------------------
1) mit 'run_Iterations_over1Signal'
Wird die Simulation mit einem Signal gestartet. Am Ende wird ein Plot 
ausgegeben auf dem der Verlauf der Korrelationskoeffizienten �ber die in 
'f_parames.m' gew�hlten Iterationen ausgegeben wird. Dabei k�nnen alle 
Iterationen angezeigt werden, wenn 'parameters.showAllIterations = 1' 
gesetzt wird. Nur die Baseline, erste und letzte Iteration angezeigt wenn
'parameters.showAllIterations = 0' ist. Es werden zudem die 1. Iteration 
und letzte Iteration der verarbeiten Signale als .wav-File abgespeichert, so 
dass sie auralisiert werden k�nnen. (Simulationsergebnisse) 

2) mit 'run_Iterations_overXSignals'
wird der arithmetische Mittelwert der Korralationskoeffizienten �ber 
X Signale, sowie deren Standardabweichung berechnet. Es werden vier 
Diagramme ausgegeben. In Figure 1 sind Baseline, erste Iteration und 
zehnte Iteration geplottet. Die anderen Diagramme stellen jeweils einen 
Korrelationsverlauf mit den dazugeh�rigen Standardabweichungen pro Frequenz-
kanal dar. Die verwendeten Signale m�ssen die gleichen Dateinamen aufweisen
und dieser fortlaufend nummeriert sein (z.B. Sprache1, Sprache2, Sprache3, 
etc.). Die gew�nschte Dateibezeichnung kann im Skript
'run_Iterations_overXSignals' ver�ndert werden. 

---------------------------------------------------------------------------
Modell und Algorithmus
---------------------------------------------------------------------------
'f_process.m' ist der Kern der Implementierung und setzt den Algorithmus des 
Modells um. Hier werden die Filterbankobjekte der NH- und 
HI-Analysefilterbanken sowie der NH- und HI-Modelle erzeugt und es findet 
die iterative Berechnung des Korrelationskoeffizienten zwischen den 
Einh�llenden statt.    

---------------------------------------------------------------------------
PATH Umgebung
---------------------------------------------------------------------------
Der Ordner 'CodeDemo' und alle seine Unterordner sind im Matlab Pfad 
einzubinden.

---------------------------------------------------------------------------
Parameter
---------------------------------------------------------------------------
Simulationsparameter werden in der Datei 'f_params.m' den entsprechenden 
Simulationsanforderungen angepasst. N�heres zu den Parametern ist 
in 'f_params' kommentiert.

---------------------------------------------------------------------------
Signals (Sinust�ne, Tonkomplexe)
---------------------------------------------------------------------------
Sollte mit Sinust�nen oder Tonkomplexen gearbeitet werden, k�nnen diese in 
den entsprechenden Zeilen der Datei 'Signals.m' konfiguriert werden.
Alle Eintr�ge in 'Signals.m', die mit Sprache, St�rger�usch und Musik 
zusammenh�ngen sollen NICHT ver�ndert werden.

---------------------------------------------------------------------------
Audio Daten
---------------------------------------------------------------------------
Der Unterordner AudioData enth�lt die Audiodateien f�r die Simulation.
Soundfiles (einkanalig bzw. mono) k�nnen hier abgelegt werden.

---------------------------------------------------------------------------
Formatierung von Diagrammen
---------------------------------------------------------------------------
Der Unterordner 'GraphicsFormatting' enth�lt kleine Funktionen mit denen 
die Diagramme am Ende einer Simulation aufbereitet werden. Diese Funktionen 
werden von den Simulationsskripten aufgerufen.

---------------------------------------------------------------------------
Sound_IO
---------------------------------------------------------------------------
Der Unterordner Sound_IO enth�lt kleine Funktionen mit denen die Diagramme 
am Ende einer Simulation aufbereitet werden. Diese Funktionen werden von 
den Simulationsskripten aufgerufen.
