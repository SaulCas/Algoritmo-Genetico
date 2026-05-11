import random
import copy
import time
import csv

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    print("Falta instalar matplotlib. Ejecuta: python -m pip install matplotlib")
    raise


# Matriz BLOSUM62 incluida para que el proyecto no dependa de la libreria blosum.
AMINOACIDOS = list("ARNDCQEGHILKMFPSTWYV")
BLOSUM62_VALORES = [
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],
    [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],
    [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],
    [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4],
]
BLOSUM62 = {
    a: {b: BLOSUM62_VALORES[i][j] for j, b in enumerate(AMINOACIDOS)}
    for i, a in enumerate(AMINOACIDOS)
}


def get_sequences():
    seq1 = "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDLYDDDDKDRWGKLVVLGAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQV"
    seq2 = "MKTLLVAAAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQKELQKQLGQKAKEL"
    seq3 = "MAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKALCVFAIN"
    return [list(seq1), list(seq2), list(seq3)]


def crear_individuo():
    return get_sequences()


def copiar_individuo(individuo):
    return [fila[:] for fila in individuo]


def crear_poblacion_inicial(n=10):
    individuo_base = crear_individuo()
    return [copiar_individuo(individuo_base) for _ in range(n)]


def igualar_longitud_secuencias(individuo, gap='-'):
    max_len = max(len(fila) for fila in individuo)
    return [fila + [gap] * (max_len - len(fila)) for fila in individuo]


def evaluar_individuo_blosum62(individuo, contador):
    contador["NFE"] += 1
    individuo = igualar_longitud_secuencias(individuo)
    score = 0
    n_seqs = len(individuo)
    seq_len = len(individuo[0])

    for col in range(seq_len):
        for i in range(n_seqs):
            for j in range(i + 1, n_seqs):
                a = individuo[i][col]
                b = individuo[j][col]
                if a == '-' or b == '-':
                    score -= 4
                else:
                    score += BLOSUM62[a][b]
    return score


def validar_poblacion_sin_gaps(poblacion, originales):
    for individuo in poblacion:
        for seq, seq_orig in zip(individuo, originales):
            seq_sin_gaps = [a for a in seq if a != '-']
            seq_orig_sin_gaps = [a for a in seq_orig if a != '-']
            if seq_sin_gaps != seq_orig_sin_gaps:
                return False
    return True


def obtener_best(scores, poblacion):
    idx_mejor = scores.index(max(scores))
    return copiar_individuo(poblacion[idx_mejor]), scores[idx_mejor]


def seleccionar_mejores(poblacion, scores, cantidad):
    indices = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    nueva_poblacion = [copiar_individuo(poblacion[i]) for i in indices[:cantidad]]
    nuevos_scores = [scores[i] for i in indices[:cantidad]]
    return nueva_poblacion, nuevos_scores


def mutar_insertando_gaps(individuo, n_gaps=1, probabilidad=1.0):
    nuevo_individuo = []
    for secuencia in individuo:
        sec = secuencia[:]
        if random.random() < probabilidad:
            posiciones = set()
            for _ in range(n_gaps):
                pos = random.randint(0, len(sec))
                while pos in posiciones:
                    pos = random.randint(0, len(sec))
                posiciones.add(pos)
                sec.insert(pos, '-')
        nuevo_individuo.append(sec)
    return nuevo_individuo


def mutar_moviendo_gaps(individuo, movimientos=1, probabilidad=0.8):
    nuevo_individuo = []
    for secuencia in individuo:
        sec = secuencia[:]
        if random.random() < probabilidad:
            for _ in range(movimientos):
                posiciones_gap = [i for i, valor in enumerate(sec) if valor == '-']
                if not posiciones_gap:
                    continue
                pos_quitar = random.choice(posiciones_gap)
                sec.pop(pos_quitar)
                pos_insertar = random.randint(0, len(sec))
                sec.insert(pos_insertar, '-')
        nuevo_individuo.append(sec)
    return nuevo_individuo


def mutar_eliminando_gaps(individuo, n_gaps=1, probabilidad=0.2):
    nuevo_individuo = []
    for secuencia in individuo:
        sec = secuencia[:]
        if random.random() < probabilidad:
            for _ in range(n_gaps):
                posiciones_gap = [i for i, valor in enumerate(sec) if valor == '-']
                if posiciones_gap:
                    sec.pop(random.choice(posiciones_gap))
        nuevo_individuo.append(sec)
    return nuevo_individuo


def cruzar_individuos_original(ind1, ind2):
    hijo1 = []
    hijo2 = []
    for seq1, seq2 in zip(ind1, ind2):
        aa_indices = [i for i, a in enumerate(seq1) if a != '-']
        if len(aa_indices) < 6:
            hijo1.append(seq1[:])
            hijo2.append(seq2[:])
            continue

        intentos = 0
        while True:
            p1, p2 = sorted(random.sample(aa_indices, 2))
            if p2 - p1 >= 5 or intentos > 10:
                break
            intentos += 1

        def cruza(seqA, seqB):
            aaA = [a for a in seqA if a != '-']
            aaB = [a for a in seqB if a != '-']
            nueva = aaA[:p1] + aaB[p1:p2] + aaA[p2:]
            resultado = []
            idx = 0
            for a in seqA:
                if a == '-':
                    resultado.append('-')
                else:
                    resultado.append(nueva[idx])
                    idx += 1
            return resultado

        hijo1.append(cruza(seq1, seq2))
        hijo2.append(cruza(seq2, seq1))

    hijo1 = mutar_insertando_gaps(hijo1, n_gaps=1, probabilidad=0.8)
    hijo2 = mutar_insertando_gaps(hijo2, n_gaps=1, probabilidad=0.8)
    return hijo1, hijo2


def cruzar_poblacion_original(poblacion):
    nueva_poblacion = []
    n = len(poblacion)
    indices = list(range(n))
    random.shuffle(indices)
    parejas = [(indices[i], indices[i + 1]) for i in range(0, n - 1, 2)]
    if n % 2 == 1:
        parejas.append((indices[-1], indices[0]))

    for idx1, idx2 in parejas:
        padre1 = poblacion[idx1]
        padre2 = poblacion[idx2]
        hijo1, hijo2 = cruzar_individuos_original(padre1, padre2)
        nueva_poblacion.append(copiar_individuo(padre1))
        nueva_poblacion.append(copiar_individuo(padre2))
        nueva_poblacion.append(hijo1)
        nueva_poblacion.append(hijo2)

    return nueva_poblacion[:2 * n]


def obtener_posiciones_gap(secuencia):
    return [i for i, valor in enumerate(secuencia) if valor == '-']


def construir_con_gaps(aminoacidos, posiciones_gap):
    posiciones_gap = sorted(set(pos for pos in posiciones_gap if pos >= 0))

    while True:
        total = len(aminoacidos) + len(posiciones_gap)
        posiciones_filtradas = [pos for pos in posiciones_gap if pos < total]
        if len(posiciones_filtradas) == len(posiciones_gap):
            break
        posiciones_gap = posiciones_filtradas

    posiciones_gap = set(posiciones_gap)
    resultado = []
    idx = 0

    for posicion in range(total):
        if posicion in posiciones_gap:
            resultado.append('-')
        else:
            resultado.append(aminoacidos[idx])
            idx += 1

    return resultado


def cruzar_individuos_mejorado(ind1, ind2):
    hijo1 = []
    hijo2 = []

    for seq1, seq2 in zip(ind1, ind2):
        aa = [a for a in seq1 if a != '-']
        gaps1 = obtener_posiciones_gap(seq1)
        gaps2 = obtener_posiciones_gap(seq2)
        max_len = max(len(seq1), len(seq2))

        if max_len <= 2:
            hijo1.append(seq1[:])
            hijo2.append(seq2[:])
            continue

        corte = random.randint(1, max_len - 1)

        gaps_hijo1 = [g for g in gaps1 if g < corte] + [g for g in gaps2 if g >= corte]
        gaps_hijo2 = [g for g in gaps2 if g < corte] + [g for g in gaps1 if g >= corte]

        hijo1.append(construir_con_gaps(aa, gaps_hijo1))
        hijo2.append(construir_con_gaps(aa, gaps_hijo2))

    return hijo1, hijo2


def seleccion_torneo(poblacion, scores, k=3):
    participantes = random.sample(range(len(poblacion)), k)
    mejor = max(participantes, key=lambda i: scores[i])
    return copiar_individuo(poblacion[mejor])


def busqueda_local(individuo, contador, intentos=4):
    mejor = igualar_longitud_secuencias(copiar_individuo(individuo))
    mejor_score = evaluar_individuo_blosum62(mejor, contador)

    for _ in range(intentos):
        candidato = copiar_individuo(mejor)
        candidato = mutar_moviendo_gaps(candidato, movimientos=1, probabilidad=0.9)
        candidato = mutar_insertando_gaps(candidato, n_gaps=1, probabilidad=0.25)
        candidato = mutar_eliminando_gaps(candidato, n_gaps=1, probabilidad=0.20)
        candidato = igualar_longitud_secuencias(candidato)
        score_candidato = evaluar_individuo_blosum62(candidato, contador)

        if score_candidato > mejor_score:
            mejor = candidato
            mejor_score = score_candidato

    return mejor, mejor_score


def ejecutar_original(generaciones=100, tam_poblacion=10, semilla=42):
    random.seed(semilla)
    contador = {"NFE": 0}
    inicio = time.time()
    historial = []

    poblacion = crear_poblacion_inicial(tam_poblacion)
    poblacion = [mutar_insertando_gaps(ind, n_gaps=1, probabilidad=1.0) for ind in poblacion]
    poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]
    scores = [evaluar_individuo_blosum62(ind, contador) for ind in poblacion]
    poblacion, scores = seleccionar_mejores(poblacion, scores, max(2, tam_poblacion // 2))

    very_best = None
    fitness_very_best = None

    for generacion in range(1, generaciones + 1):
        poblacion = cruzar_poblacion_original(poblacion)
        poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]
        scores = [evaluar_individuo_blosum62(ind, contador) for ind in poblacion]
        poblacion, scores = seleccionar_mejores(poblacion, scores, max(2, tam_poblacion // 2))
        best, fitness_best = obtener_best(scores, poblacion)

        if very_best is None or fitness_best > fitness_very_best:
            very_best = best
            fitness_very_best = fitness_best

        historial.append({
            "algoritmo": "Original",
            "generacion": generacion,
            "fitness": fitness_very_best,
            "NFE": contador["NFE"],
            "tiempo": time.time() - inicio
        })

    validacion = validar_poblacion_sin_gaps([very_best], get_sequences())
    return historial, very_best, fitness_very_best, contador["NFE"], time.time() - inicio, validacion


def ejecutar_mejorado(generaciones=100, tam_poblacion=24, semilla=42):
    random.seed(semilla)
    contador = {"NFE": 0}
    inicio = time.time()
    historial = []
    elite = 4

    poblacion = []
    for _ in range(tam_poblacion):
        individuo = crear_individuo()
        individuo = mutar_insertando_gaps(individuo, n_gaps=random.randint(1, 3), probabilidad=1.0)
        individuo = mutar_moviendo_gaps(individuo, movimientos=random.randint(0, 2), probabilidad=0.8)
        poblacion.append(igualar_longitud_secuencias(individuo))

    scores = [evaluar_individuo_blosum62(ind, contador) for ind in poblacion]
    very_best, fitness_very_best = obtener_best(scores, poblacion)

    for generacion in range(1, generaciones + 1):
        poblacion, scores = seleccionar_mejores(poblacion, scores, tam_poblacion)
        nueva_poblacion = [copiar_individuo(ind) for ind in poblacion[:elite]]

        while len(nueva_poblacion) < tam_poblacion:
            padre1 = seleccion_torneo(poblacion, scores, k=3)
            padre2 = seleccion_torneo(poblacion, scores, k=3)
            hijo1, hijo2 = cruzar_individuos_mejorado(padre1, padre2)

            hijo1 = mutar_moviendo_gaps(hijo1, movimientos=2, probabilidad=0.85)
            hijo1 = mutar_insertando_gaps(hijo1, n_gaps=1, probabilidad=0.35)
            hijo1 = mutar_eliminando_gaps(hijo1, n_gaps=1, probabilidad=0.18)
            hijo1, _ = busqueda_local(hijo1, contador, intentos=3)
            nueva_poblacion.append(hijo1)

            if len(nueva_poblacion) < tam_poblacion:
                hijo2 = mutar_moviendo_gaps(hijo2, movimientos=2, probabilidad=0.85)
                hijo2 = mutar_insertando_gaps(hijo2, n_gaps=1, probabilidad=0.35)
                hijo2 = mutar_eliminando_gaps(hijo2, n_gaps=1, probabilidad=0.18)
                hijo2, _ = busqueda_local(hijo2, contador, intentos=3)
                nueva_poblacion.append(hijo2)

        poblacion = [igualar_longitud_secuencias(ind) for ind in nueva_poblacion]
        scores = [evaluar_individuo_blosum62(ind, contador) for ind in poblacion]
        best, fitness_best = obtener_best(scores, poblacion)

        if fitness_best > fitness_very_best:
            very_best = best
            fitness_very_best = fitness_best

        historial.append({
            "algoritmo": "Mejorado",
            "generacion": generacion,
            "fitness": fitness_very_best,
            "NFE": contador["NFE"],
            "tiempo": time.time() - inicio
        })

    validacion = validar_poblacion_sin_gaps([very_best], get_sequences())
    return historial, very_best, fitness_very_best, contador["NFE"], time.time() - inicio, validacion


def guardar_csv(historial_completo, nombre_archivo="resultados_original_vs_mejorado.csv"):
    with open(nombre_archivo, "w", newline="", encoding="utf-8") as archivo:
        columnas = ["algoritmo", "generacion", "fitness", "NFE", "tiempo"]
        escritor = csv.DictWriter(archivo, fieldnames=columnas)
        escritor.writeheader()
        escritor.writerows(historial_completo)


def graficar(historial_original, historial_mejorado):
    gen_o = [dato["generacion"] for dato in historial_original]
    fit_o = [dato["fitness"] for dato in historial_original]
    nfe_o = [dato["NFE"] for dato in historial_original]
    time_o = [dato["tiempo"] for dato in historial_original]

    gen_m = [dato["generacion"] for dato in historial_mejorado]
    fit_m = [dato["fitness"] for dato in historial_mejorado]
    nfe_m = [dato["NFE"] for dato in historial_mejorado]
    time_m = [dato["tiempo"] for dato in historial_mejorado]

    plt.figure(figsize=(10, 6))
    plt.plot(gen_o, fit_o, label="Algoritmo original")
    plt.plot(gen_m, fit_m, label="Algoritmo mejorado")
    plt.title("Comparación de fitness: original vs mejorado")
    plt.xlabel("Generación")
    plt.ylabel("Mejor fitness acumulado")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("comparacion_original_vs_mejorado.png", dpi=150)
    plt.close()

    plt.figure(figsize=(10, 6))
    plt.plot(gen_o, nfe_o, label="Algoritmo original")
    plt.plot(gen_m, nfe_m, label="Algoritmo mejorado")
    plt.title("Comparación de NFE: original vs mejorado")
    plt.xlabel("Generación")
    plt.ylabel("Número de evaluaciones")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("comparacion_nfe_original_vs_mejorado.png", dpi=150)
    plt.close()

    plt.figure(figsize=(10, 6))
    plt.plot(gen_o, time_o, label="Algoritmo original")
    plt.plot(gen_m, time_m, label="Algoritmo mejorado")
    plt.title("Comparación de tiempo: original vs mejorado")
    plt.xlabel("Generación")
    plt.ylabel("Tiempo acumulado en segundos")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("comparacion_tiempo_original_vs_mejorado.png", dpi=150)
    plt.close()


def guardar_resumen(fitness_original, fitness_mejorado, nfe_original, nfe_mejorado, tiempo_original, tiempo_mejorado, val_original, val_mejorado):
    mejora = fitness_mejorado - fitness_original
    porcentaje = (mejora / abs(fitness_original)) * 100 if fitness_original != 0 else 0
    texto = f"""RESUMEN DEL PROYECTO FINAL

Resultado del algoritmo original:
Fitness final: {fitness_original}
NFE final: {nfe_original}
Tiempo final: {tiempo_original:.4f} segundos
Validacion de integridad: {val_original}

Resultado del algoritmo mejorado:
Fitness final: {fitness_mejorado}
NFE final: {nfe_mejorado}
Tiempo final: {tiempo_mejorado:.4f} segundos
Validacion de integridad: {val_mejorado}

Mejora obtenida:
Diferencia de fitness: {mejora}
Porcentaje aproximado de mejora: {porcentaje:.2f}%

Mejoras aplicadas:
1. Poblacion inicial mas diversa.
2. Elitismo para conservar los mejores individuos.
3. Seleccion por torneo para elegir padres con mejor fitness.
4. Cruza enfocada en intercambiar posiciones de gaps.
5. Mutacion mejorada que puede insertar, mover y eliminar gaps.
6. Busqueda local para intentar mejorar cada hijo antes de pasarlo a la siguiente generacion.
7. Validacion de integridad para comprobar que las secuencias siguen siendo las mismas al quitar los gaps.
"""
    with open("resumen_proyecto_final.txt", "w", encoding="utf-8") as archivo:
        archivo.write(texto)


def main():
    generaciones = 100
    semilla = 42

    print("Ejecutando algoritmo original...")
    historial_original, best_original, fitness_original, nfe_original, tiempo_original, val_original = ejecutar_original(
        generaciones=generaciones,
        tam_poblacion=10,
        semilla=semilla
    )

    print("Ejecutando algoritmo mejorado...")
    historial_mejorado, best_mejorado, fitness_mejorado, nfe_mejorado, tiempo_mejorado, val_mejorado = ejecutar_mejorado(
        generaciones=generaciones,
        tam_poblacion=24,
        semilla=semilla
    )

    historial_completo = historial_original + historial_mejorado
    guardar_csv(historial_completo)
    graficar(historial_original, historial_mejorado)
    guardar_resumen(
        fitness_original,
        fitness_mejorado,
        nfe_original,
        nfe_mejorado,
        tiempo_original,
        tiempo_mejorado,
        val_original,
        val_mejorado
    )

    print("\nRESULTADOS FINALES")
    print("Algoritmo original -> fitness:", fitness_original, "NFE:", nfe_original, "tiempo:", round(tiempo_original, 4), "validacion:", val_original)
    print("Algoritmo mejorado -> fitness:", fitness_mejorado, "NFE:", nfe_mejorado, "tiempo:", round(tiempo_mejorado, 4), "validacion:", val_mejorado)
    print("Mejora de fitness:", fitness_mejorado - fitness_original)

    print("\nArchivos generados:")
    print("comparacion_original_vs_mejorado.png")
    print("comparacion_nfe_original_vs_mejorado.png")
    print("comparacion_tiempo_original_vs_mejorado.png")
    print("resultados_original_vs_mejorado.csv")
    print("resumen_proyecto_final.txt")


if __name__ == "__main__":
    main()
