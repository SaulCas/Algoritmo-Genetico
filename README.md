# Proyecto Final: Versión mejorada de Algoritmo Genético

Este proyecto compara un algoritmo genético original contra una versión mejorada.

## Objetivo

Crear una versión propia y mejorada del algoritmo genético para intentar aumentar el fitness, manteniendo la validación de integridad de las secuencias.

## Archivos principales

- `AG10_proyecto_final_mejorado.py`: código principal del proyecto.
- `comparacion_original_vs_mejorado.png`: gráfica principal de fitness.
- `comparacion_nfe_original_vs_mejorado.png`: comparación de evaluaciones NFE.
- `comparacion_tiempo_original_vs_mejorado.png`: comparación de tiempo.
- `resultados_original_vs_mejorado.csv`: resultados por generación.
- `resumen_proyecto_final.txt`: resumen automático de resultados.

## Cómo ejecutar

Instalar matplotlib si hace falta:

```bash
python -m pip install matplotlib
```

Ejecutar:

```bash
python AG10_proyecto_final_mejorado.py
```

## Mejoras aplicadas

1. Población inicial más diversa.
2. Elitismo para conservar los mejores individuos.
3. Selección por torneo para elegir padres con mejor fitness.
4. Cruza enfocada en intercambiar posiciones de gaps.
5. Mutación mejorada que puede insertar, mover y eliminar gaps.
6. Búsqueda local para mejorar hijos antes de pasar a la siguiente generación.
7. Conservación de la validación de integridad.

## Validación de integridad

La validación revisa que, al quitar los gaps (`-`), las secuencias resultantes sigan siendo iguales a las secuencias originales.

Si el programa imprime:

```text
validacion: True
```
