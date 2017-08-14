# This Makefile is not normally needed when installing Muscato using
# go get.  If you modify the source code in your local go/src
# directory, you can use make to recompile and install all the
# components of Muscato.

all: muscato muscato_screen muscato_confirm muscato_prep_reads muscato_prep_targets muscato_window_reads

.PHONY: muscato muscato_screen muscato_confirm muscato_prep_reads muscato_prep_targets muscato_window_reads

muscato:
	go install ./muscato

muscato_screen:
	go install ./muscato_screen

muscato_confirm:
	go install ./muscato_confirm

muscato_prep_reads:
	go install ./muscato_prep_reads

muscato_prep_targets:
	go install ./muscato_prep_targets

muscato_window_reads:
	go install ./muscato_window_reads
