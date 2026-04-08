from __future__ import annotations
import os, re, shutil
from typing import Optional, Tuple

BEGIN = "# BEGIN STRESSVIZ LOVE"
END   = "# END STRESSVIZ LOVE"

def _fmt_complex(z: complex) -> str:
    sign = "-" if z.imag < 0 else "+"
    return f"{z.real:.6g} {sign} {abs(z.imag):.6g}i"

def _parse_complex(s: str) -> complex:
    s = s.strip().replace("−","-").replace("i","j").replace("I","j")
    if s.startswith(("+","-")) and "j" in s and not re.search(r"[0-9]\s*[+\-]", s[1:]):
        s = "0+" + s
    return complex(s)

def read_love_from_sat(path: str) -> Optional[Tuple[complex, complex, complex]]:
    if not os.path.isfile(path):
        return None
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            txt = f.read()
    except Exception:
        return None

    m = re.search(rf"(?ms)^{re.escape(BEGIN)}\s*(.*?)^\s*{re.escape(END)}\s*$", txt)
    if not m:
        return None
    body = m.group(1)

    def grab(key_options):
        for k in key_options:
            mm = re.search(rf"(?mi)^\s*#?\s*{k}\s*[:=]\s*(.+?)\s*$", body)
            if mm:
                return _parse_complex(mm.group(1))
        return None

    h2 = grab(("h2","h"))
    k2 = grab(("k2","k"))
    l2 = grab(("l2","l"))
    if None in (h2, k2, l2):
        return None
    return h2, k2, l2

def write_love_to_sat(path: str, h2: complex, k2: complex, l2: complex, make_backup: bool = True) -> bool:
    if not os.path.isfile(path):
        return False
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        txt = f.read()

    block = (
        f"{BEGIN}\n"
        f"# h2: {_fmt_complex(h2)}\n"
        f"# k2: {_fmt_complex(k2)}\n"
        f"# l2: {_fmt_complex(l2)}\n"
        f"{END}\n"
    )

    if re.search(rf"(?ms)^{re.escape(BEGIN)}.*?{re.escape(END)}\s*$", txt):
        new_txt = re.sub(rf"(?ms)^{re.escape(BEGIN)}.*?{re.escape(END)}\s*$", block.strip()+"\n", txt)
    else:
        sep = "" if txt.endswith("\n") else "\n"
        new_txt = txt + sep + block

    if make_backup:
        try:
            shutil.copyfile(path, path + ".bak")
        except Exception:
            pass

    with open(path, "w", encoding="utf-8") as f:
        f.write(new_txt)
    return True
