// Copyright 2017, Kerby Shedden and the Muscato contributors.

package utils

func CountDinuc(seq []byte, wk []int) int {

	for i, _ := range wk {
		wk[i] = 0
	}

	var last int
	var n int
	for i, x := range seq {

		var v int
		switch x {
		case 'A':
			v = 0
		case 'T':
			v = 1
		case 'G':
			v = 2
		case 'C':
			v = 3
		default:
			v = 4
		}

		if i > 0 {
			k := 5*last + v
			if wk[k] == 0 {
				n++
			}
			wk[k]++
		}
		last = v
	}

	return n
}
