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
	filePath  string
	schema    *arrow.Schema
	writer    *ipc.FileWriter
	chunkSize int
	mu        sync.Mutex
}

func NewArrowWriter(filePath string, fieldNames []string, fieldTypes []arrow.DataType, chunkSize int) (*ArrowWriter, error) {
	schema := makeSchema(fieldNames, fieldTypes)

	file, err := os.Create(filePath)
	if err != nil {
		return nil, err
	}

	writer, err := ipc.NewFileWriter(file, ipc.WithSchema(schema))
	if err != nil {
		return nil, err
	}

	return &ArrowWriter{
		filePath:  filePath,
		schema:    schema,
		writer:    writer,
		chunkSize: chunkSize,
	}, nil
}

func (aw *ArrowWriter) writeChunk(record arrow.Record) error {
	aw.mu.Lock()
	defer aw.mu.Unlock()

	if err := aw.writer.Write(record); err != nil {
		return err
	}

	return nil
}

func (aw *ArrowWriter) Close() error {
	aw.mu.Lock()
	defer aw.mu.Unlock()

	return aw.writer.Close()
}

type ArrowRowBuilder struct {
	builders       []array.Builder
	appendFuncs    []func(array.Builder, interface{}) error
	pool           *memory.GoAllocator
	arrowWriter    *ArrowWriter
	numRowsInChunk int
}

func NewArrowRowBuilder(aw *ArrowWriter) (*ArrowRowBuilder, error) {
	pool := memory.NewGoAllocator()
	builders, appendFuncs, err := makeBuilders(aw.schema, pool)

	if err != nil {
		return nil, err
	}

	return &ArrowRowBuilder{
		builders:       builders,
		appendFuncs:    appendFuncs,
		pool:           pool,
		arrowWriter:    aw,
		numRowsInChunk: 0,
	}, nil
}

func (arb *ArrowRowBuilder) WriteRow(row []interface{}) error {
	if len(row) != len(arb.builders) {
		return fmt.Errorf("mismatch in number of fields: expected %d, got %d", len(arb.builders), len(row))
	}

	for i, val := range row {
		if err := arb.appendFuncs[i](arb.builders[i], val); err != nil {
			return fmt.Errorf("error appending to column %d: %v", i, err)
		}
	}

	arb.numRowsInChunk++

	if arb.numRowsInChunk == arb.arrowWriter.chunkSize {
		arb.writeChunk()
		arb.numRowsInChunk = 0
	}

	return nil
}

func (arb *ArrowRowBuilder) writeChunk() error {
	var cols []arrow.Array
	for _, b := range arb.builders {
		// NewArray creates a new array from the builder and resets the builder
		// See: https://github.com/apache/arrow/blob/maint-14.0.2/go/arrow/array/builder.go#L84
		cols = append(cols, b.NewArray())
	}

	record := array.NewRecord(arb.arrowWriter.schema, cols, int64(arb.numRowsInChunk))
	defer record.Release()

	return nil
}

func (arb *ArrowRowBuilder) Release() error {
	if arb.numRowsInChunk > 0 {
		err := arb.writeChunk()
		arb.numRowsInChunk = 0

		if err != nil {
			return err
		}
	}

	for _, builder := range arb.builders {
		builder.Release()
	}

	return nil
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

func makeSchema(fieldNames []string, fieldTypes []arrow.DataType) *arrow.Schema {
	fields := make([]arrow.Field, len(fieldTypes))
	for i, dataType := range fieldTypes {
		fields[i] = arrow.Field{Name: fieldNames[i], Type: dataType}
	}

	schema := arrow.NewSchema(fields, nil)

	schema.Fields()

	return schema
}

func makeBuilders(schema *arrow.Schema, pool *memory.GoAllocator) ([]array.Builder, []func(array.Builder, interface{}) error, error) {
	builders := make([]array.Builder, schema.NumFields())
	appendFuncs := make([]func(array.Builder, interface{}) error, schema.NumFields())

	for i, field := range schema.Fields() {
		switch field.Type {
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
			return nil, nil, fmt.Errorf("unsupported data type: %s", field.Type)
		}
	}
	return builders, appendFuncs, nil
}
