#!/usr/bin/env bash
# Download the three Fedkiw fire papers into ./papers
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p papers
curl -L -o papers/2001_nguyen_fedkiw_kang_flame_discontinuities.pdf \
  https://physbam.stanford.edu/~fedkiw/papers/cam2000-19.pdf
curl -L -o papers/2002_nguyen_fedkiw_jensen_fire.pdf \
  https://physbam.stanford.edu/papers/stanford2002-02.pdf
curl -L -o papers/2007_hong_shinar_fedkiw_wrinkled_flames.pdf \
  https://physbam.stanford.edu/papers/stanford2007-02.pdf
echo "done → papers/"
ls -la papers/
