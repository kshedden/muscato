all: muscato muscato_screen muscato_confirm muscato_prep_reads muscato_prep_targets muscato_window_reads

muscato: muscato/muscato.go
	go install ./muscato

muscato_screen: muscato_screen/muscato_screen.go
	go install ./muscato_screen

muscato_confirm: muscato_confirm/muscato_confirm.go
	go install ./muscato confirm

muscato_prep_reads: muscato_prep_reads/muscato_prep_reads.go
	go install ./muscato_prep_reads

muscato_prep_targets: muscato_prep_targets/muscato_prep_targets.go
	go install ./muscato_prep_targets

muscato_window_reads: muscato_window_reads/muscato_window_reads.go
	go install ./muscato_window_reads
