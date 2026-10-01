# PET MCP Geant4

Simulation Geant4 minimale d'un détecteur PET formé d'empilements de plaques à
microcanaux (MCP) placés de part et d'autre de la source. La géométrie (nombre
de plaques de chaque côté, dimensions, canaux) se configure depuis une macro. Chaque événement émet deux photons
gamma de 511 keV dos à dos. Le programme enregistre les électrons qui atteignent
un canal MCP et les coïncidences entre les deux côtés.

## Architecture

```text
main.cc
 ├─ DetectorConstruction     empilements de MCP en verre plombé et leurs canaux
 │   └─ DetectorMessenger     commandes macro /mcp/...
 ├─ FTFP_BERT                physique de référence
 │   └─ Livermore + Auger    physique électromagnétique basse énergie
 └─ ActionInitialization
     ├─ PrimaryGeneratorAction   deux gammas dos à dos
     │   └─ PrimaryGeneratorMessenger  commandes macro /source/...
     ├─ SteppingAction           détection à l'entrée des canaux
     ├─ EventAction              multiplicité et coïncidence par événement
     └─ RunAction                résumé du run et fichier ROOT
```

La simulation s'arrête à l'entrée d'un électron dans un canal : elle ne modélise
pas l'avalanche électronique interne au MCP. Par défaut, le transport du gamma
s'arrête à sa sortie de la plaque (`/mcp/killGammaOnExit`). Ces deux conditions
font partie du modèle courant.

## Modèle physique

- monde et canaux : vide `G4_Galactic` ;
- plaques : verre `GlassLead60Perc`, densité 4,93 g/cm³, fractions massiques
  Pb/Si/O de 60/16/24 % ;
- rayon d'une plaque : 25 mm ; épaisseur : 1,5 mm ;
- canaux : diamètre 25 µm, pas hexagonal 31 µm, inclinaison ±8° ;
- cut de production dans les MCP : 1 µm ;
- source : deux gammas de 511 keV à l'origine, opposés, avec direction uniforme
  en angle solide dans un cône de demi-angle 75° autour de `+z`, entièrement
  reconfigurable par macro ;
- physique : `FTFP_BERT`, physique EM Livermore et émission Auger activée.

## Configuration de la géométrie

Toutes les dimensions sont pilotables par macro via les commandes `/mcp/...`.
`macros/geometry.mac` décrit la géométrie de référence et est exécutée par les
autres macros juste avant `/run/initialize`. Une macro appelante peut
surcharger n'importe laquelle de ses lignes après l'avoir exécutée :

| commande | effet | défaut |
| --- | --- | --- |
| `/mcp/stacksPlus <n>` | nombre de plaques empilées côté `+z` (`0` supprime le côté) | `1` |
| `/mcp/stacksMinus <n>` | nombre de plaques empilées côté `-z` | `1` |
| `/mcp/stacks <n>` | raccourci : même nombre des deux côtés | — |
| `/mcp/firstPlateZ <d> <u>` | distance origine → centre de la première plaque | `1 cm` |
| `/mcp/stackGap <d> <u>` | vide entre deux plaques consécutives | `100 um` |
| `/mcp/plateThickness <d> <u>` | épaisseur d'une plaque | `1500 um` |
| `/mcp/plateRadius <d> <u>` | rayon d'une plaque | `25 mm` |
| `/mcp/channelDiameter <d> <u>` | diamètre d'un canal | `25 um` |
| `/mcp/channelPitch <d> <u>` | pas hexagonal | `31 um` |
| `/mcp/channelAngle <d> <u>` | inclinaison des canaux (signe suivant le côté) | `8 deg` |
| `/mcp/chevron true\|false` | alterne le signe de l'inclinaison d'une plaque à l'autre | `false` |
| `/mcp/edgeMarginFraction <f>` | bord sans canaux, en fraction du rayon | `0.2` |
| `/mcp/worldSize <d> <u>` | côté du monde cubique, agrandi d'office s'il est trop petit | `10 cm` |
| `/mcp/productionCut <d> <u>` | cut de production dans la région MCP | `1 um` |
| `/mcp/killGammaOnExit true\|false` | arrête le gamma dès qu'il quitte sa plaque | `true` |
| `/mcp/print` | affiche la configuration courante | — |

Exemple : trois plaques chevron par côté, gamma laissé libre de traverser une
plaque pour atteindre la suivante.

```text
/mcp/stacks 3
/mcp/chevron true
/mcp/stackGap 200 um
/mcp/killGammaOnExit false
```

Les commandes sont acceptées avant `/run/initialize` et entre deux runs : dans
ce second cas la géométrie est reconstruite automatiquement avant le `beamOn`
suivant.

`/mcp/killGammaOnExit` conserve par défaut la règle historique : un gamma est
tué dès qu'il sort de sa plaque. Avec plus d'une plaque par côté il faut la
désactiver, sinon le gamma n'atteint jamais la plaque suivante.

## Configuration de la source

Même principe pour la source, avec les commandes `/source/...` regroupées dans
`macros/source.mac` :

| commande | effet | défaut |
| --- | --- | --- |
| `/source/particle <nom>` | particule émise | `gamma` |
| `/source/energy <e> <u>` | énergie de chaque particule | `511 keV` |
| `/source/backToBack true\|false` | émet une seconde particule opposée depuis le même vertex | `true` |
| `/source/perEvent <n>` | nombre d'annihilations (de paires si `backToBack`) par événement | `1` |
| `/source/position <x> <y> <z> <u>` | centre du point d'émission | `0 0 0 mm` |
| `/source/positionSigma <d> <u>` | étalement gaussien du vertex sur chaque axe (`0` = ponctuelle) | `0 mm` |
| `/source/coneAxis <ux> <uy> <uz>` | axe du cône d'émission | `0 0 1` |
| `/source/coneHalfAngle <a> <u>` | demi-angle du cône ; `180 deg` = isotrope | `75 deg` |
| `/source/print` | affiche la configuration courante | — |

Le tirage reste uniforme en angle solide dans le cône. Comme pour `/mcp/...`,
ces commandes sont acceptées avant `/run/initialize` et entre deux runs.
`/source/coneHalfAngle 180 deg` donne une source isotrope.

## Compilation

Prérequis : CMake 3.16+, un compilateur C++20, Geant4 avec `ui_all` et `vis_all`,
et ROOT avec `Core`, `RIO` et `Tree`.

```bash
cmake -S . -B build
cmake --build build --parallel
```

## Exécution

CMake copie les macros dans `build/macros/`. Exécutez l'application depuis le
dossier de build afin d'y regrouper les sorties générées :

```bash
cd build
./PET_MCP macros/run.mac
```

`run.mac` utilise une graine fixe et lance 1 000 événements. Le résultat est
écrit dans `build/mcp_output.root`.

Pour la visualisation :

```bash
cd build
./PET_MCP
```

Sans argument, l'application ouvre la session interactive et charge
`macros/vis.mac`. Les millions de canaux paramétrés rendent l'affichage très
lourd ; pour vérifier simplement l'agencement, réduisez la plaque avant
`/run/initialize`, par exemple `/mcp/plateRadius 0.5 mm` et
`/mcp/channelPitch 60 um`.

Les macros disponibles :

| macro | rôle |
| --- | --- |
| `geometry.mac` | géométrie de référence, chargée par toutes les autres |
| `source.mac` | source de référence, chargée par toutes les autres |
| `run.mac` | 1 000 événements, graine fixe |
| `test.mac` | run court de 10 événements |
| `run_photons.mac` | étude longue de 500 000 événements |
| `vis.mac` | visualisation de la géométrie de référence |
| `scan.mac` | balayage du nombre de plaques par côté |
| `scan_point.mac` | un point du balayage, appelé par `scan.mac` |

Pour lancer un run avec une autre configuration, remplacez la ligne
`/control/execute geometry.mac` ou `/control/execute source.mac` de la macro de
run, ou ajoutez les commandes `/mcp/...` et `/source/...` voulues juste après.

## Sortie ROOT et analyse

`ElectronChannelHitTree` contient une ligne lors de l'entrée d'un électron dans
un canal : événement, track, côté, rang de la plaque dans l'empilement
(`stackIndex`, compté depuis l'origine), numéro du canal, énergie, temps,
position, direction et processus créateur.

`ScanTree`, présent uniquement dans le fichier de balayage, contient une ligne
par run (voir « Balayage de paramètres »).

`EventSummaryTree` contient une ligne par événement : nombres de hits côté `+z`
et `-z`, ainsi que l'indicateur de coïncidence. Les événements sans hit sont
conservés pour définir correctement l'efficacité.

Pour lancer l'étude d'un million de photons :

```bash
cd build
./PET_MCP macros/run_photons.mac
cd ..
root -l 'analysis/AnalyzePhoton.C("build/mcp_output.root")'
```

La macro ROOT affiche le taux de photons ayant produit au moins un électron qui
atteint un canal, puis trace le spectre énergétique de ces électrons. Elle ne
crée aucun fichier ROOT supplémentaire.

## Balayage de paramètres

`macros/scan.mac` balaie le nombre de plaques MCP par côté (symétrique) avec
`/control/foreach`, et `macros/scan_point.mac` décrit un point du balayage.
Chaque run ajoute une ligne à l'arbre `ScanTree` d'un fichier cumulatif unique,
créé par `/output/scanFile` :

```bash
cd build
./PET_MCP macros/scan.mac
cd ..
root -l 'analysis/AnalyzeStackScan.C("build/scan_stacks.root")'
```

`AnalyzeStackScan.C` affiche un tableau par point (événements, électrons
atteignant un canal, répartition `+z` / `-z`, coïncidences) et trace les
électrons collectés en fonction du nombre de plaques dans `scan_stacks.png`.

Colonnes de `ScanTree` : `stacksPlus`, `stacksMinus`, `events`,
`electronHits`, `electronHitsPlus`, `electronHitsMinus`, `eventsWithHitPlus`,
`eventsWithHitMinus`, `coincidences`.

Les commandes de sortie :

| commande | effet | défaut |
| --- | --- | --- |
| `/output/file <nom>` | fichier détaillé, recréé à chaque `beamOn` | `mcp_output.root` |
| `/output/scanFile <nom>` | démarre un balayage : crée le fichier et y ajoute une ligne par run | — |

Pour balayer un autre paramètre, remplacez la ligne `/mcp/stacks {stacks}` de
`scan_point.mac` par la commande voulue : la mécanique `foreach` et le fichier
cumulatif ne changent pas. `scan.mac` fixe `/mcp/killGammaOnExit false` pour
tous les points, afin que seule la profondeur de l'empilement varie.

Le balayage à la géométrie nominale est coûteux (~1,5 M de canaux par plaque,
soit ~12 M de copies au point à 4 plaques par côté). `scan.mac` contient une
ligne `/mcp/plateRadius` commentée pour une mise au point rapide.

## Arborescence

```text
.
├── CMakeLists.txt
├── README.md
├── main.cc              point d'entrée
├── include/             en-têtes
├── src/                 implémentations
├── macros/              macros Geant4
└── analysis/            macros ROOT d'analyse
```

## Limites

L'application utilise un RunManager série. Elle ne simule ni le gain MCP, ni
l'électronique de lecture, ni la reconstruction PET. La géométrie compte plus de
trois millions de copies paramétrées au total ; sa construction domine le temps
d'un run court.

Aucune licence n'est actuellement définie. Ce point doit être réglé avant toute
publication ou redistribution.
