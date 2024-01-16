package arrow

import (
	"fmt"
	"os"
	"sync"

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
	appendFuncs    []func(array.Builder, interface{}) error
	mu             sync.Mutex
}

func NewArrowWriter(filePath string, fieldNames []string, fieldTypes []arrow.DataType, chunkSize int) (*ArrowWriter, error) {
	pool := memory.NewGoAllocator()
	fields := make([]arrow.Field, len(fieldTypes))
	builders := make([]array.Builder, len(fieldTypes))
	appendFuncs := make([]func(array.Builder, interface{}) error, len(fieldTypes))

	for i, dataType := range fieldTypes {
		fields[i] = arrow.Field{Name: fieldNames[i], Type: dataType}
		switch dataType {
		case arrow.PrimitiveTypes.Uint16:
			builders[i] = array.NewUint16Builder(pool)
			appendFuncs[i] = appendUint16
		case arrow.PrimitiveTypes.Uint32:
			builders[i] = array.NewUint32Builder(pool)
			appendFuncs[i] = appendUint32
		case arrow.PrimitiveTypes.Uint64:
			builders[i] = array.NewUint64Builder(pool)
			appendFuncs[i] = appendUint64
		case arrow.PrimitiveTypes.Int16:
			builders[i] = array.NewInt16Builder(pool)
			appendFuncs[i] = appendInt16
		case arrow.PrimitiveTypes.Int32:
			builders[i] = array.NewInt32Builder(pool)
			appendFuncs[i] = appendInt32
		case arrow.PrimitiveTypes.Int64:
			builders[i] = array.NewInt64Builder(pool)
			appendFuncs[i] = appendInt64
		case arrow.PrimitiveTypes.Float32:
			builders[i] = array.NewFloat32Builder(pool)
			appendFuncs[i] = appendFloat32
		case arrow.PrimitiveTypes.Float64:
			builders[i] = array.NewFloat64Builder(pool)
			appendFuncs[i] = appendFloat64
		case arrow.BinaryTypes.String:
			builders[i] = array.NewStringBuilder(pool)
			appendFuncs[i] = appendString
		case arrow.FixedWidthTypes.Boolean:
			builders[i] = array.NewBooleanBuilder(pool)
			appendFuncs[i] = appendBool
		default:
			return nil, fmt.Errorf("unsupported data type: %s", dataType)
		}
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

	return &ArrowWriter{
		filePath:       filePath,
		schema:         schema,
		writer:         writer,
		builders:       builders,
		pool:           pool,
		chunkSize:      chunkSize,
		numRowsInChunk: 0,
		appendFuncs:    appendFuncs,
	}, nil
}

func (aw *ArrowWriter) Write(row []interface{}) error {
	aw.mu.Lock()
	defer aw.mu.Unlock()

	if len(row) != len(aw.builders) {
		return fmt.Errorf("mismatch in number of fields: expected %d, got %d", len(aw.builders), len(row))
	}

	for i, val := range row {
		if err := aw.appendFuncs[i](aw.builders[i], val); err != nil {
			return fmt.Errorf("error appending to column %d: %v", i, err)
		}
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
	aw.mu.Lock()
	defer aw.mu.Unlock()

	// Write any remaining data
	if aw.numRowsInChunk > 0 {
		if err := aw.writeChunk(); err != nil {
			return err
		}
	}
	return aw.writer.Close()
}

func appendUint16(builder array.Builder, val interface{}) error {
	if v, ok := val.(uint16); ok {
		builder.(*array.Uint16Builder).Append(v)
		return nil
	}
	return fmt.Errorf("type mismatch, expected uint16")
}

func appendUint32(builder array.Builder, val interface{}) error {
	if v, ok := val.(uint32); ok {
		builder.(*array.Uint32Builder).Append(v)
		return nil
	}
	return fmt.Errorf("type mismatch, expected uint32")
}

func appendUint64(builder array.Builder, val interface{}) error {
	if v, ok := val.(uint64); ok {
		builder.(*array.Uint64Builder).Append(v)
		return nil
	}
	return fmt.Errorf("type mismatch, expected uint64")
}

func appendInt16(builder array.Builder, val interface{}) error {
	if v, ok := val.(int16); ok {
		builder.(*array.Int16Builder).Append(v)
		return nil
	}

	return fmt.Errorf("type mismatch, expected int16")
}

func appendInt32(builder array.Builder, val interface{}) error {
	if v, ok := val.(int32); ok {
		builder.(*array.Int32Builder).Append(v)
		return nil
	}

	return fmt.Errorf("type mismatch, expected int32")
}

func appendInt64(builder array.Builder, val interface{}) error {
	if v, ok := val.(int64); ok {
		builder.(*array.Int64Builder).Append(v)
		return nil
	}

	return fmt.Errorf("type mismatch, expected int64")
}

func appendFloat32(builder array.Builder, val interface{}) error {
	if v, ok := val.(float32); ok {
		builder.(*array.Float32Builder).Append(v)
		return nil
	}

	return fmt.Errorf("type mismatch, expected float32")
}

func appendFloat64(builder array.Builder, val interface{}) error {
	if v, ok := val.(float64); ok {
		builder.(*array.Float64Builder).Append(v)
		return nil
	}

	return fmt.Errorf("type mismatch, expected float64")
}

func appendString(builder array.Builder, val interface{}) error {
	if v, ok := val.(string); ok {
		builder.(*array.StringBuilder).Append(v)
		return nil
	}

	return fmt.Errorf("type mismatch, expected string")
}

func appendBool(builder array.Builder, val interface{}) error {
	if v, ok := val.(bool); ok {
		builder.(*array.BooleanBuilder).Append(v)
		return nil
	}

	return fmt.Errorf("type mismatch, expected bool")
}
