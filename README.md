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

Il raggio di ogni particella viene stimato in funzione della massa, senza superare un massimo calcolato a partire dalla distanza dei due corpi più vicini moltiplicata per un fattore `0.1`. Tutti gli altri raggi vengono calcolati in proporzione alla massa più grande ed elevati alla ⅓: la radice cubica serve per avvicinare i raggi tra loro — ad esempio, un corpo avente massa pari alla metà di quella massima avrà un raggio circa `0.8` volte il raggio massimo. Questo evita che due corpi si sovrappongano già nelle condizioni iniziali e collidano al primo step, e che un corpo enormemente più massivo degli altri copra tutta la finestra.

---

## Come eseguire il programma

Dopo aver configurato e buildato il progetto, dalla cartella del progetto:

```bash
cmake -S . -B build -G"Ninja Multi-Config"
cmake --build build --config Debug
./build/Debug/nbody input/nome_file.txt
```

dove `input/nome_file.txt` è il path del file di configurazione passato come argomento al `main()`.

---

## Librerie esterne

L'unica libreria esterna da installare è SFML:

```bash
sudo apt install libsfml-dev
```

---

## Output e implementazione grafica

L'output avviene su due interfacce: terminale e finestra grafica.

Sul **terminale** viene stampato ogni secondo:
- Numero di step computati
- Energia meccanica
- Energia cinetica
- Energia potenziale
- Quantità di moto
- Momento angolare

Queste quantità si conservano durante la simulazione (tranne le energie in caso di collisioni).

Tramite **SFML** è implementata la visualizzazione grafica delle orbite. La simulazione gira a 60 FPS. Premendo la barra spaziatrice è possibile fermare e far ripartire la simulazione.

---

## Gestione dello zoom

- Se `zoom > 0.` nel file di input, lo zoom è fisso a quel valore
- Se `zoom = 0.`, il programma calcola automaticamente lo zoom per ogni frame:
  1. Trova il corpo più lontano dall'origine (`max_d`)
  2. Calcola lo zoom ideale come `350.f / max_d` (350 px è poco meno della metà della larghezza della finestra): maggiore è la distanza dall'origine, minore è lo zoom, così che tutti i corpi rimangano visibili
  3. Applica uno smoothing: invece di scattare al valore ideale, ogni frame aggiunge il 5% della differenza tra zoom corrente e zoom ideale, mantenendo una continuità di visualizzazione

In questo modo è possibile visualizzare correttamente sistemi con masse di ordini di grandezza molto diversi senza dover impostare uno zoom predefinito.

---

## Strategia e creazione dei test

I test sono divisi in cinque fasi:

1. **Verifica algebrica** — operazioni su `Vec` tra vettori e con scalari, eccezioni come massa negativa
2. **Verifica analitica delle funzioni fisiche** — `computePosition()`, `computeVelocity()`, `totalAcceleration()`, energie e momenti confrontati con valori calcolati analiticamente
3. **Rendimento dinamico e collisioni** — coerenza dei dati prima e dopo il Velocity Verlet, corretta gestione delle collisioni
4. **Conservazione delle quantità fisiche** — energia meccanica, momento angolare e quantità di moto dopo N step
5. **Verifica dell'input** — `read_from()`: file inesistente, file vuoto, commenti, riga `@`, dati fisicamente impossibili

---

## Interpretazione dei risultati

Le simulazioni risultano coerenti con i dati fisici attesi. Come esempio tangibile, il file `solar_system_reduced.txt` simula il sistema solare ridotto: le orbite dei pianeti sono stabili e le quantità fisiche (energia meccanica, quantità di moto, momento angolare) si conservano entro oscillazioni accettabili, confermando la correttezza dell'integratore simplettico Velocity Verlet.

> **Nota:** l'unica imprecisione riscontrata è che in seguito a un urto non sempre viene conservato il momento angolare totale del sistema — fenomeno fisicamente errato, probabilmente dovuto al fatto che non è stato considerato il contributo dello spin della nuova massa generata dalla fusione.

---

## Uso di sistemi di Intelligenza Artificiale Generativa

Sono stati utilizzati sistemi di intelligenza artificiale come strumento di documentazione e per generare esempi di file di input, e più raramente come supporto per agevolare la comprensione dei messaggi del compilatore.
