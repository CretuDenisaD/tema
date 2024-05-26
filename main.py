import tkinter as tk
from tkinter import messagebox
import cmath
import math
import matplotlib.pyplot as plt


def complexToString(z):
    x = round(z.real, 3)
    y = round(z.imag, 3)
    if y == 0:
        return f"{x:g}"
    if x == 0:
        return f"{y:g}j"
    return f"{x:g}{y:+g}j"


def rezolvareEcuatie(a, b, c, d):
    def Radix3(w):
        rho, theta = cmath.polar(w)
        return [cmath.rect(math.pow(rho, 1 / 3), (theta + 2 * k * math.pi) / 3) for k in range(3)]

    assert a != 0, "a must be non-zero"
    p = c / a - b * b / (3.0 * a * a)
    q = 2.0 * b ** 3 / (27.0 * a ** 3) - b * c / (3.0 * a * a) + d / a
    Delta = q * q + 4.0 * p ** 3 / 27.0
    w = (-q + cmath.sqrt(Delta)) / 2
    results = []
    solutions = []
    for k, u in enumerate(Radix3(w)):
        z = u - p / (3 * u)
        z -= b / (3.0 * a)
        solutions.append(z)
        results.append(f"z{k} = {complexToString(z)}")
    return results, solutions


def VerficareDeltaMinus(p, q):
    delta = q * q + 4.0 * p ** 3 / 27.0
    assert delta < 0, "Nu suntem in cazul Delta < 0"
    theta = math.acos(3 * q * math.sqrt(-3 / p) / (2 * p))
    rho = math.sqrt(-p / 3)
    results = []
    solutions = []
    results.append(f"Ecuatia z**3{p:+}z{q:+}=0")
    results.append("are solutiile:")
    for k in range(3):
        zk = 2 * rho * math.cos((theta + 2 * k * math.pi) / 3)
        solutions.append(zk)
        results.append(f"z{k} = {zk:.3g}")
    return results, solutions


def solve_ecuatie():
    try:
        a = complex(a_entry.get())
        b = complex(b_entry.get())
        c = complex(c_entry.get())
        d = complex(d_entry.get())
        results, solutions = rezolvareEcuatie(a, b, c, d)
        result_label.config(text="\n".join(results))
        plot_solutions(solutions, "Solutions of Cubic Equation")
    except Exception as e:
        messagebox.showerror("Error", str(e))


def solve_verificare():
    try:
        p = float(p_entry.get())
        q = float(q_entry.get())
        results, solutions = VerficareDeltaMinus(p, q)
        verificare_result_label.config(text="\n".join(results))
        plot_real_solutions(solutions, "Real Solutions for Delta < 0")
    except Exception as e:
        messagebox.showerror("Error", str(e))


def plot_solutions(solutions, title):
    plt.figure()
    for sol in solutions:
        plt.plot(sol.real, sol.imag, 'ro')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.title(title)
    plt.xlabel('Real Part')
    plt.ylabel('Imaginary Part')
    plt.show()


def plot_real_solutions(solutions, title):
    plt.figure()
    for sol in solutions:
        plt.plot(sol, 0, 'ro')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.title(title)
    plt.xlabel('Real Part')
    plt.ylabel('Imaginary Part')
    plt.show()


root = tk.Tk()
root.title("Rezolvare ecuatii grad 3")
root.geometry("500x700")

title_label = tk.Label(root, text="Rezolvare ecuatii de grad 3")
title_label.pack()

frame_ecuatie = tk.Frame(root)
frame_ecuatie.pack(pady=10)

valori_label = tk.Label(frame_ecuatie, text="Valori pentru ecuatia de grad 3:")
valori_label.grid(row=0, column=0, padx=5, pady=5, columnspan=2)

a_label = tk.Label(frame_ecuatie, text="a:")
a_label.grid(row=1, column=0, padx=5, pady=5)
a_entry = tk.Entry(frame_ecuatie)
a_entry.grid(row=1, column=1, padx=5, pady=5)

b_label = tk.Label(frame_ecuatie, text="b:")
b_label.grid(row=2, column=0, padx=5, pady=5)
b_entry = tk.Entry(frame_ecuatie)
b_entry.grid(row=2, column=1, padx=5, pady=5)

c_label = tk.Label(frame_ecuatie, text="c:")
c_label.grid(row=3, column=0, padx=5, pady=5)
c_entry = tk.Entry(frame_ecuatie)
c_entry.grid(row=3, column=1, padx=5, pady=5)

d_label = tk.Label(frame_ecuatie, text="d:")
d_label.grid(row=4, column=0, padx=5, pady=5)
d_entry = tk.Entry(frame_ecuatie)
d_entry.grid(row=4, column=1, padx=5, pady=5)

solve_button = tk.Button(frame_ecuatie, text="Rezolva", command=solve_ecuatie)
solve_button.grid(row=5, column=0, columnspan=2, pady=10)

result_label = tk.Label(root, text="")
result_label.pack(pady=20)

frame_verificare = tk.Frame(root)
frame_verificare.pack(pady=10)

verificare_label = tk.Label(frame_verificare, text="Verificare Delta Minus:")
verificare_label.grid(row=0, column=0, padx=5, pady=5, columnspan=2)

p_label = tk.Label(frame_verificare, text="p:")
p_label.grid(row=1, column=0, padx=5, pady=5)
p_entry = tk.Entry(frame_verificare)
p_entry.grid(row=1, column=1, padx=5, pady=5)

q_label = tk.Label(frame_verificare, text="q:")
q_label.grid(row=2, column=0, padx=5, pady=5)
q_entry = tk.Entry(frame_verificare)
q_entry.grid(row=2, column=1, padx=5, pady=5)

verificare_button = tk.Button(frame_verificare, text="Verifica", command=solve_verificare)
verificare_button.grid(row=3, column=0, columnspan=2, pady=10)

verificare_result_label = tk.Label(root, text="")
verificare_result_label.pack(pady=20)

root.mainloop()
