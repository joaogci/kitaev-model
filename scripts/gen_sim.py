output_file = "sims"

def frange(start, stop, step):
  while start <= stop:
    yield start
    start += step


N = 182 # 182 364 448
alpha = 1.0

T_vals = [round(j, 2) for j in list(frange(0.01, 2.05, 0.05))]

def run(T):
  return True


n_sim = 0
with open(output_file, "w") as f:
  for T_, T in enumerate(T_vals):
    run_str = "N"
    if run(T):
      run_str = "Y"
      n_sim += 1

    line = f"{run_str}  {N}  {alpha}  {round(1.0/T, 2)}   50  10 \n"
    f.write(line)

  f.write("stop\n")
  f.write("Run N_R alpha_R  Beta_R  NSW_R NB_R\n")
  f.write("Name=(Run N_R alpha_R  Beta_R  NSW_R NB_R)\n")

print(f"#simulations: {n_sim}")
