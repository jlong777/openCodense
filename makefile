CC=gcc -std=gnu99 -Wall
CFLAGS=-O3
DEBUG=-g
OMP=-fopenmp
PROF=#-pg

codense: CreateMicroarrayDataSetGraph.o CreateSummaryGraph.o                    \
         ODES_master.o GenerateEdgeSupportMatrix.o codense.h                    \
         GenerateSecondOrderGraph.o ReadMicroarrayData.o Vertex2Edge.o          \
         util.o ODES_slave.o CodenseMI.o codense.c
	$(CC) $(CFLAGS) $(DEBUG) $(OMP) $(PROF) codense.c CodenseMI.o           \
	CreateMicroarrayDataSetGraph.o CreateSummaryGraph.o ODES_master.o       \
	Vertex2Edge.o GenerateEdgeSupportMatrix.o GenerateSecondOrderGraph.o    \
	ReadMicroarrayData.o ODES_slave.o util.o -o codense -lm -lpthread

CodenseMI.o: CodenseMI.c codense.h
	$(CC) -c $(CFLAGS) $(DEBUG) $(PROF) CodenseMI.c
CreateMicroarrayDataSetGraph.o: CreateMicroarrayDataSetGraph.c codense.h
	$(CC) -c $(CFLAGS) $(DEBUG) $(OMP) $(PROF) CreateMicroarrayDataSetGraph.c
CreateSummaryGraph.o: CreateSummaryGraph.c codense.h
	$(CC) -c $(CFLAGS) $(DEBUG) $(PROF) CreateSummaryGraph.c
ODES_master.o: ODES_master.c codense.h
	$(CC) -c $(CFLAGS) $(DEBUG) $(PROF) ODES_master.c
ODES_slave.o: ODES_slave.c codense.h
	$(CC) -c           $(DEBUG) $(PROF) ODES_slave.c
GenerateEdgeSupportMatrix.o: GenerateEdgeSupportMatrix.c codense.h
	$(CC) -c $(CFLAGS) $(DEBUG) $(PROF) GenerateEdgeSupportMatrix.c
GenerateSecondOrderGraph.o: GenerateSecondOrderGraph.c codense.h
	$(CC) -c $(CFLAGS) $(DEBUG) $(PROF) GenerateSecondOrderGraph.c
ReadMicroarrayData.o: ReadMicroarrayData.c codense.h
	$(CC) -c $(CFLAGS) $(DEBUG) $(PROF) ReadMicroarrayData.c
Vertex2Edge.o: Vertex2Edge.c codense.h
	$(CC) -c $(CFLAGS) $(DEBUG) $(PROF) Vertex2Edge.c
util.o: util.c codense.h
	$(CC) -c $(CFLAGS) $(DEBUG) $(PROF) util.c

clean: 
	rm *o

distclean:
	if [ -e codense ];    then rm codense;    fi
	if [ -e parameters ]; then rm parameters; fi
	rm *o
        
        
# Copyright (C) 2015 James Long
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
