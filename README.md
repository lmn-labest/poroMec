# poroMec

Repositorio do modelo poromecanio 3D.

Index
- [poroMec](#poromec)
- [Pre](#pre)
  - [Compilar o metis](#compilar-o-metis)
  - [Compilar o pre](#compilar-o-pre)
- [Rodando o pre](#rodando-o-pre)

# Pre

## Compilar o metis

Descompactar e compilar o metis5

```console
tar -xvzf metis-5.1.0.tar.gz
cd metis-5.1.0
make
```

É necessario colocar a lib do metis na pasta `lib` do projeto do pre

```console
lib
└── libmetis.a
```

## Compilar o pre

Fazer uma copia do `Makefile` base

```console
cp Makefile_base Makefile
```

As versões novas do metis precisam dessa lib extra `libGKlib`. As versões mais antigas do metis podem ser encontrados [aqui](http://glaros.dtc.umn.edu/gkhome/metis/metis/download).

Para compilar basta:

```console
make
```

# Rodando o pre

Criar um arquivo `pre.dat` com o conteudo

```
input  mesh.dat
output     part
div          12
method non-overllaping
partVtk     yes
partMeshVtk yes
partMeshMef yes
meshLoads   yes
vtkBin       no
memory     1000
end
```

Para rodar basta

```
prepar pre.dar
```
