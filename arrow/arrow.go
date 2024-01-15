package arrow

import (
	"fmt"
	"os"

	"github.com/apache/arrow/go/v14/arrow"
	"github.com/apache/arrow/go/v14/arrow/array"
	"github.com/apache/arrow/go/v14/arrow/ipc"
	"github.com/apache/arrow/go/v14/arrow/memory"
)

type ArrowWriter struct {
	filePath       string
	schema         *arrow.Schema
	writer         *ipc.FileWriter
	builders       []array.Builder
	pool           *memory.GoAllocator
	chunkSize      int
	numRowsInChunk int
}

func NewArrowWriter(filePath string, fieldNames []string, chunkSize int) (*ArrowWriter, error) {
	pool := memory.NewGoAllocator()
	fields := make([]arrow.Field, len(fieldNames))

	for i, name := range fieldNames {
		fields[i] = arrow.Field{Name: name, Type: arrow.PrimitiveTypes.Float64}
	}

	schema := arrow.NewSchema(fields, nil)
	file, err := os.Create(filePath)
	if err != nil {
		return nil, err
	}

	writer, err := ipc.NewFileWriter(file, ipc.WithSchema(schema))
	if err != nil {
		return nil, err
	}

	builders := make([]array.Builder, len(fields))
	for i := range fields {
		builders[i] = array.NewFloat64Builder(pool)
	}

	return &ArrowWriter{
		filePath:       filePath,
		schema:         schema,
		writer:         writer,
		builders:       builders,
		pool:           pool,
		chunkSize:      chunkSize,
		numRowsInChunk: 0,
	}, nil
}

func (aw *ArrowWriter) Write(row []float64) error {
	if len(row) != len(aw.builders) {
		return fmt.Errorf("mismatch in number of fields: expected %d, got %d", len(aw.builders), len(row))
	}

	for i, val := range row {
		aw.builders[i].(*array.Float64Builder).Append(val)
	}

	aw.numRowsInChunk++

	if aw.numRowsInChunk == aw.chunkSize {
		if err := aw.writeChunk(); err != nil {
			return err
		}
	}

	return nil
}

func (aw *ArrowWriter) writeChunk() error {
	var cols []arrow.Array
	for _, b := range aw.builders {
		// NewArray creates a new array from the builder and resets the builder
		// See: https://github.com/apache/arrow/blob/maint-14.0.2/go/arrow/array/builder.go#L84
		cols = append(cols, b.NewArray())
	}

	record := array.NewRecord(aw.schema, cols, int64(aw.numRowsInChunk))
	defer record.Release()

	if err := aw.writer.Write(record); err != nil {
		return err
	}

	// Reset for the next chunk
	aw.numRowsInChunk = 0

	return nil
}

func (aw *ArrowWriter) Close() error {
	// Write any remaining data
	if aw.numRowsInChunk > 0 {
		if err := aw.writeChunk(); err != nil {
			return err
		}
	}
	return aw.writer.Close()
}
