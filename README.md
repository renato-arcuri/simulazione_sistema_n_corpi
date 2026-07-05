# Simulazione Sistema a N-Corpi
**Arcuri Renato - Brollo Nicolò - Collini Emanuele**  
*03/07/2026*

---

## Implementazione degli N-corpi

I corpi nella simulazione sono stati implementati tramite una classe `Particle`, avente come membri privati la massa, la posizione, la velocità, il raggio e l'accelerazione. Questi ultimi due vengono calcolati a partire dagli altri, mentre gli altri vengono inizializzati in input. Il raggio è calcolato in funzione della massa, mentre l'accelerazione viene ricalcolata a ogni step di simulazione in base all'interazione gravitazionale. Le grandezze fisiche vettoriali (posizione, velocità e accelerazione) sono di tipo `Vec`, una struct contenente 2 `double`. Il sistema è quindi un `std::vector<Particle>`.

---

## Implementazione Velocity Verlet

Ogni step temporale è eseguito dalla funzione `step()`.  
Per ogni oggetto `Particle` contenuto nel vettore `system` vengono eseguiti in ordine i seguenti passaggi:

1. Calcolo e modifica delle posizioni: `computePosition()`
2. Calcolo delle accelerazioni: `setAcceleration()` attraverso `totalAcceleration()`, che calcola l'accelerazione che subisce la particella i-esima rispetto alle altre presenti
3. Calcolo e modifica delle velocità: `computeVelocity()`

Infine viene verificata l'eventuale presenza di collisioni tramite `manageCollisions()`: se la distanza tra i centri di due particelle è minore o uguale alla somma dei rispettivi raggi, si elimina la seconda particella mentre le proprietà della prima vengono aggiornate in accordo alle leggi fisiche che regolano un urto totalmente anelastico.

---

## Input

L'input avviene tramite la lettura di un file di configurazione. In esso sono immessi per ogni particella le condizioni di partenza: la massa, la posizione iniziale e la velocità iniziale. Vengono inoltre forniti i parametri di configurazione della simulazione: la costante di gravitazione universale `G`, lo zoom, lo step temporale `dt` e la `speed`.

- Se `zoom` è impostato a `0.` il programma lo calcolerà automaticamente; se è un valore diverso, sarà fisso durante tutta l'esecuzione
- `speed` è il numero di volte che viene iterato `step()` per ogni frame (~17 ms a 60 FPS); maggiore è il valore, maggiore è la velocità di simulazione

Ogni riga del file (che non inizi per `@` o `#`) rappresenta una particella. La riga dei parametri di configurazione deve avere `@` come primo carattere, con valori in ordine: `G zoom dt speed`. Questa riga è opzionale; se assente vengono usati valori predefiniti. Ogni valore è gestito come `double` (tranne `speed`, che è un `int`), è quindi preferibile esprimerlo con il punto anche se il valore non ha cifre significative dopo la virgola (es. `0.`). Le righe che iniziano con `#` sono commenti e vengono ignorate.

Il formato di ogni riga particella è: `mass rx ry vx vy`

Esempio di file di input:

```text
# figure-8 di Chenciner-Montgomery, G=1
# G   zoom  dt    speed
@ 1.0 200.  0.001 17
# mass rx           ry          vx           vy
1.  -0.97000436  0.24308753  0.4662036850  0.4323657300
1.   0.97000436 -0.24308753  0.4662036850  0.4323657300
1.   0.          0.         -0.93240737   -0.86473146
```

Nella cartella `./input` del progetto sono presenti diversi esempi:

- `figure8.txt`
- `solar_system_reduced.txt`
- `solar_system.txt`
- `scattering.txt`
- `spinning_system.txt`
- `collision.txt`
- `caotic_system.txt`

---

## Calcolo del raggio a partire dalla massa

Il raggio di ogni particella viene stimato in funzione della massa, senza superare un massimo calcolato a partire dalla distanza dei due corpi più vicini moltiplicata per un fattore `0.1`. Tutti gli altri raggi vengono calcolati in proporzione alla massa più grande ed elevati alla ⅓: la radice
