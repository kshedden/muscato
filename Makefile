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
