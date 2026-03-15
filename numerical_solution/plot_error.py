from pathlib import Path

import matplotlib.pyplot as plt


DATA_FILE = Path("data/error_history.dat")
OUTPUT_FILE = Path("data/error_history.png")


def main() -> None:
    if not DATA_FILE.exists():
        raise FileNotFoundError(f"Не найден файл с историей ошибки: {DATA_FILE}")

    iterations = []
    errors = []

    with DATA_FILE.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            iterations.append(int(parts[0]))
            errors.append(float(parts[1]))

    if not iterations:
        raise ValueError("В файле нет данных для построения графика.")

    plt.figure(figsize=(8, 5))
    plt.plot(iterations, errors, marker="o", linewidth=1.5)
    plt.yscale("log")
    plt.xlabel("Iteration")
    plt.ylabel("Error")
    plt.title("Error vs Iteration")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(OUTPUT_FILE, dpi=200)
    plt.show()


if __name__ == "__main__":
    main()
