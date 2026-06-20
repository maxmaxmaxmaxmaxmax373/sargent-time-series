# Maintainer's Guide — Editing & Publishing This Book

This book is written in **MyST Markdown** (plus a couple of Jupyter notebooks) and
published as a website with **Jupyter Book** + **GitHub Pages**. This guide takes
you from a fresh download all the way to "my change is live online."

> **Mental model.** MyST markdown = the recipe. `jupyter-book build` = the cooking.
> GitHub Actions = a robot chef that cooks automatically on every push.
> GitHub Pages = the restaurant that serves the dish at a `github.io` address.

---

## 0. What's in this folder

```
.
├── _config.yml          # book settings (title, math macros, repo URL)
├── _toc.yml             # table of contents — which chapters appear, in what order
├── requirements.txt     # the tools needed to build (jupyter-book, numpy, scipy, matplotlib)
├── intro.md             # the landing page
├── chapters/            # all 40 chapters: *.md (markdown) and *.ipynb (notebooks)
├── code/                # Python scripts that generate the figures
├── figures/             # the pre-rendered figure PNGs that chapters embed
├── data/                # cached data (e.g. FRED downloads)
└── .github/workflows/
    └── deploy.yml       # the "robot chef" — auto-builds & publishes on every push
```

---

## 1. Install the prerequisites (one time)

You need three things on your computer:

| Tool | Check it's installed | If missing |
|------|----------------------|------------|
| **Git** | `git --version` | https://git-scm.com/downloads |
| **Python 3.11+** | `python3 --version` | https://www.python.org/downloads |
| **A GitHub account** | — | https://github.com/join |

> The cloud build uses Python **3.11**; matching it locally avoids surprises.

---

## 2. Get the project onto your computer

**If you already have the folders** (you downloaded them), open a terminal *inside*
this folder (the one containing `_config.yml`).

**Or clone a fresh copy** from GitHub:

```bash
git clone https://github.com/maxmaxmaxmaxmaxmax373/sargent-time-series.git
cd sargent-time-series
```

To publish later you must be able to push to a GitHub repo. Two options:

- **A — You have write access** to `maxmaxmaxmaxmaxmax373/sargent-time-series`:
  nothing more to do; skip to Step 3.
- **B — You want your *own* online copy** (recommended if you're not the owner):
  see [Appendix: Publish your own copy](#appendix-publish-your-own-copy) at the end,
  then come back.

---

## 3. Set up the build environment (one time per machine)

Create an isolated Python environment and install the tools listed in
`requirements.txt`:

```bash
# from inside the book folder (where _config.yml lives)
python3 -m venv .venv

# activate it:
source .venv/bin/activate        # macOS / Linux
# .venv\Scripts\activate         # Windows (PowerShell)

pip install -r requirements.txt
```

You'll know it worked if `jupyter-book --version` prints a version number.

> **Every new terminal session**, re-activate the environment with the
> `source .venv/bin/activate` line before building. (You only `pip install` once.)

---

## 4. Build & preview the book locally

Before publishing anything, build it on your own machine and look at it:

```bash
jupyter-book build .
```

This writes a website into `_build/html/`. Open the front page in your browser:

```bash
open _build/html/index.html          # macOS
# xdg-open _build/html/index.html    # Linux
# start _build/html/index.html       # Windows
```

Click around. **This local copy is private to you** — nobody else sees it until you
push. If the build prints `build succeeded`, you're good. (Warnings are usually
harmless; a red `ERROR` means something needs fixing — see
[Troubleshooting](#troubleshooting).)

> Tip: if a change doesn't show up, force a clean rebuild with
> `jupyter-book build . --all`.

---

## 5. Make your changes

### Edit existing text
Open any file in `chapters/` (e.g. `chapters/06_spectrum.md`) in a text editor,
change the prose or math, save. Math uses LaTeX between `$...$` (inline) or
`$$...$$` (display).

> ⚠️ **Display-math rule:** always leave a **blank line** before and after a
> `$$...$$` block. Gluing it to surrounding text makes the renderer swallow your
> words. ✅ Good:
> ```
> Some text.
>
> $$ y = x^2 $$
>
> More text.
> ```

### Add a brand-new chapter
1. Create `chapters/41_my_topic.md`.
2. Add it to `_toc.yml` under the right `chapters:` list:
   ```yaml
       - file: chapters/41_my_topic
   ```
   (no `.md` extension).

### Add or change a figure
Figures are **pre-generated PNGs** in `figures/`, produced by the scripts in
`code/`. To regenerate one after editing its script:
```bash
python code/fig3_kuznets_gain.py     # writes figures/fig3_kuznets_gain.png
```
Embed it in a chapter with:
```
```{figure} ../figures/fig3_kuznets_gain.png
:name: fig-3
:width: 90%

**Figure 3.** Your caption here.
```
```

### Edit a notebook chapter (the `*.ipynb` labs)
Notebook outputs are **saved into the file** and the website shows those saved
outputs (the build does not re-run notebooks). So if you change a notebook's code,
**re-run it and save** so the committed outputs match:
```bash
jupyter nbconvert --to notebook --execute --inplace chapters/40_nonlinear_lab.ipynb
```

After any change, **rebuild and preview again** (Step 4).

---

## 6. Publish: commit and push

Once the local preview looks right, send it to GitHub. Publishing is automatic
from there.

```bash
git add -A                                   # stage all your changes
git commit -m "Describe what you changed"    # save a labeled snapshot
git push                                      # upload to GitHub
```

That `git push` is the only action needed — **you do not build or upload the
website yourself.** The push triggers the robot.

---

## 7. What happens after you push (automatic)

The file `.github/workflows/deploy.yml` tells GitHub: *"on every push to `main`,
build the book and publish it."* Within ~1–2 minutes GitHub will:

1. spin up a clean Linux machine,
2. `pip install -r requirements.txt`,
3. run `jupyter-book build .`,
4. upload `_build/html/` to GitHub Pages.

**Watch it happen:** on your repo page click the **Actions** tab. A yellow dot =
building; a green ✔ = published. (A red ✘ means the build failed — click it to
read the log; the most common cause is a typo that builds fine for you only because
of a stale cache, so try `jupyter-book build . --all` locally to reproduce.)

---

## 8. Open the published book

There are **two different links** — this trips everyone up:

| Link | Shows |
|------|-------|
| `https://github.com/<user>/<repo>` | the **source files** (ingredients) — *not* the readable book |
| `https://<user>.github.io/<repo>/` | the **actual book** — formatted pages, math, figures |

👉 **Read the book at the `github.io` address:**

**https://maxmaxmaxmaxmaxmax373.github.io/sargent-time-series/**

The repo's right-sidebar **"About → 🔗 website"** field also links straight to it.
After a successful Actions run, refresh that page (use **Cmd/Ctrl + Shift + R** to
bypass the browser cache) to see your change.

---

## Quick reference (the everyday loop)

```bash
source .venv/bin/activate              # 1. activate environment
# ...edit files in chapters/ ...       # 2. make changes
jupyter-book build .                   # 3. build locally
open _build/html/index.html            # 4. preview
git add -A && git commit -m "msg"      # 5. commit
git push                               # 6. publish (auto-deploys)
# 7. refresh the github.io URL after ~2 min
```

---

## Troubleshooting

| Symptom | Fix |
|---------|-----|
| `jupyter-book: command not found` | You didn't activate the venv: `source .venv/bin/activate`. |
| Change not visible locally | `jupyter-book build . --all` (forces a full rebuild). |
| Change not visible online | Wait for the green ✔ in the **Actions** tab, then hard-refresh (Cmd/Ctrl+Shift+R). |
| Math text "eats" the words around it | Add blank lines around the `$$...$$` block (see Step 5). |
| Build fails only online (red ✘) | Open the Actions log; reproduce locally with `jupyter-book build . --all`. |
| Notebook shows old figures | Re-execute it: `jupyter nbconvert --to notebook --execute --inplace chapters/<name>.ipynb`, then commit. |

---

## Appendix: Publish your *own* copy

If you are not the owner of `maxmaxmaxmaxmaxmax373/sargent-time-series` but want
your own online edition:

1. **Create a new empty repo** on GitHub (e.g. `my-time-series-book`). Don't add a
   README.
2. **Point this folder at it** and push:
   ```bash
   git remote remove origin                       # drop the old remote (if present)
   git remote add origin https://github.com/<YOUR-USER>/my-time-series-book.git
   git branch -M main
   git push -u origin main
   ```
3. **Update the repo URL** in `_config.yml` so the in-page "edit/repository" buttons
   point to your repo:
   ```yaml
   repository:
     url: https://github.com/<YOUR-USER>/my-time-series-book
     branch: main
   ```
   Commit and push that change too.
4. **Turn on GitHub Pages with Actions:** on your repo go to
   **Settings → Pages → Build and deployment → Source**, and choose
   **"GitHub Actions"**. (The workflow file `.github/workflows/deploy.yml` is
   already included, so nothing else to configure.)
5. **Push anything** (or re-run the workflow from the Actions tab). After the green
   ✔, your book is live at:
   ```
   https://<YOUR-USER>.github.io/my-time-series-book/
   ```

That's it — from then on your everyday loop is exactly the
[Quick reference](#quick-reference-the-everyday-loop) above.
