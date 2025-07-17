# 🏔️ Diversité fonctionnelle des pelouses alpines – Mini-projet ECOMONT 2024

Ce dépôt rassemble les données, scripts et figures associés à un projet réalisé dans le cadre du **Master 2 Écologie des milieux de montagne (ECOMONT)**, porté par l’Université Savoie Mont Blanc et l’UGA. Le travail a été mené en septembre 2024 dans les pelouses alpines du **site d’Aravo** (près du col du Lautaret).

## 🧭 Objectif du projet

Explorer les liens entre les **traits fonctionnels des plantes** et les **conditions environnementales microlocales** dans une pelouse alpine.

---

## 🧪 Méthodologie

- 🌿 **Traits fonctionnels étudiés** : SLA, masse des graines, hauteur végétative, angle foliaire, épaisseur, teneur en azote, etc.
- 🌱 **Données floristiques** : 408 photos annotées, 21 espèces dominantes identifiées.
- 🌍 **Variables environnementales** : pH, humidité, matière organique, GDD, convexité, SWI, résistivité, déjections de moutons.
- 📊 **Analyses statistiques** :
  - Analyse en composantes principales (ACP)
  - RLQ & Fourth-corner
  - Krigeage spatial
  - Clustering des communautés (Bray-Curtis + NMDS)

---

## 📈 Résultats clés

- 📌 4 communautés de traits fonctionnels différenciées par des variables comme le **GDD**, le **pH**, le **SWI** ou l’**humidité**.
- 🌿 Certaines espèces investissent dans des **stratégies conservatrices** (feuilles épaisses, graines lourdes) selon le microclimat.
- 💧 Les zones humides et plates accueillent des espèces à **faible teneur en azote** et à **croissance lente**.
- 🌐 Confirmation partielle des résultats obtenus 20 ans plus tôt par Choler (2005).

---

## 🗃️ Structure du dépôt

```bash
├── scripts/              # Scripts d'analyse R
├── figures/              # Figures du rapport
├── rapport/              # Rapport PDF du projet
└── README.md             # Ce fichier
