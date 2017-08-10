all: muscato muscato_screen muscato_confirm muscato_prep_reads muscato_prep_targets muscato_window_reads

muscato: muscato/muscato.go
	go build muscato

muscato_screen: muscato_screen/muscato_screen.go
	go build muscato_screen

muscato_confirm: muscato_confirm/muscato_confirm.go
	go build muscato confirm

muscato_prep_reads: muscato_prep_reads/muscato_prep_reads.go
	go build muscato_prep_reads

muscato_prep_targets: muscato_prep_targets/muscato_prep_targets.go
	go build muscato_prep_targets

muscato_window_reads: muscato_window_reads/muscato_window_reads.go
	go build muscato_window_reads
