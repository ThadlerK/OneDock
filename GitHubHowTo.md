# Git Workflow Guide: Heidelberg_Team_2

Dieser Guide beschreibt unseren Standard-Prozess für die Arbeit mit Git. Bitte halte dich an diese Schritte, um Konflikte zu vermeiden.

---

## 1. Eigene Arbeit speichern (Täglicher Rhythmus)
Sobald du an deinen Dateien gearbeitet hast und den Stand sichern willst:

1.  **In den Projektordner wechseln:**
    ```bash
    cd Heidelberg_Team_2
    ```

2.  **Auf deinen persönlichen Branch wechseln:**
    *(Ersetze `dein-branch-name` mit deinem tatsächlichen Branch)*
    ```bash
    git checkout dein-branch-name
    ```

3.  **Änderungen vorbereiten (Stagen) & Speichern (Commit):**
    ```bash
    git add .
    git commit -m "Beschreibe kurz was du gemacht hast"
    ```

4.  **Hochladen (Push):**
    ```bash
    git push origin dein-branch-name
    ```

---

## 2. Deinen Branch aktuell halten (Updates vom Main holen)
Bevor du weiterarbeitest (oder bevor du einen Pull Request machst), solltest du die Änderungen der anderen in deinen Branch holen.

1.  **Wissen über Online-Änderungen holen (ohne Mischen):**
    ```bash
    git fetch origin
    ```

2.  **Den Main in deinen Branch mischen:**
    *(Stelle sicher, dass du auf deinem Branch bist)*
    ```bash
    git merge origin/main
    ```

3.  **Falls Konflikte auftreten (VS Code meldet "Merge Conflict"):**
    * Öffne die betroffenen Dateien (in VS Code rot markiert).
    * Entscheide dich bei den markierten Blöcken (Dein Code vs. Incoming Code).
    * Speichere die Datei(en).
    * Führe den Merge-Commit aus:
        ```bash
        git add .
        git commit -m "Konflikte mit Main gelöst"
        ```

---

## 3. Arbeit abschließen (Merge in Main)
Wenn dein Feature fertig ist und für alle im Team verfügbar sein soll.

### Option A: Über GitHub (Empfohlen)
1.  Gehe auf die GitHub-Seite des Repositories.
2.  Erstelle einen **Pull Request** (`dein-branch` ➡️ `main`).
3.  Klicke auf **Merge**.

### Option B: Manuell über das Terminal
*(Nur nutzen, wenn Option A nicht möglich ist)*

```bash
# 1. Zum Main wechseln
git checkout main

# 2. Main aktualisieren
git pull origin main

# 3. Deinen Branch hineinholen
git merge dein-branch-name

# 4. Ergebnis hochladen
git push origin main

```
