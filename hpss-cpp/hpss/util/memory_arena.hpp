#pragma once

#include <cassert>
#include <cstddef>
#include <span>

#include <chowdsp_fft.h>

namespace hpss
{
/**
 * Returns the next pointer with a given byte alignment,
 * or the base pointer if it is already aligned.
 */
template <typename Type, typename Integer_Type>
Type* snap_pointer_to_alignment (Type* base_pointer,
                                 Integer_Type alignment_bytes) noexcept
{
    return (Type*) ((((size_t) base_pointer) + (alignment_bytes - 1)) & ~(alignment_bytes - 1));
}

static constexpr int default_byte_alignment = 16;

/**
 * A simple memory arena, backed by a span of bytes.
 *
 * By default, the arena owns its memory. Alternatively,
 * set_memory_resource() points the arena at externally-owned memory.
 */
class Memory_Arena
{
public:
    Memory_Arena() = default;

    /** Constructs the arena with an initial allocated size, owned by this arena. */
    explicit Memory_Arena (size_t size_in_bytes) { reset (size_in_bytes); }

    Memory_Arena (const Memory_Arena&) = delete;
    Memory_Arena& operator= (const Memory_Arena&) = delete;

    Memory_Arena (Memory_Arena&& other) noexcept { *this = std::move (other); }
    Memory_Arena& operator= (Memory_Arena&& other) noexcept
    {
        free_owned_memory();
        raw_data = other.raw_data;
        bytes_used = other.bytes_used;
        owns_memory = other.owns_memory;
        other.raw_data = {};
        other.bytes_used = 0;
        other.owns_memory = false;
        return *this;
    }

    ~Memory_Arena() { free_owned_memory(); }

    /** (Re-)allocates the internal buffer with a given number of bytes, owned by this arena. */
    void reset (size_t new_size_bytes)
    {
        free_owned_memory();
        auto* memory = static_cast<std::byte*> (chowdsp::fft::aligned_malloc (new_size_bytes));
        raw_data = { memory, new_size_bytes };
        owns_memory = true;
        bytes_used = 0;
    }

    /**
     * Points this arena at externally-owned memory.
     * The caller is responsible for managing the lifetime of this memory.
     */
    void set_memory_resource (std::span<std::byte> external_memory) noexcept
    {
        free_owned_memory();
        raw_data = external_memory;
        owns_memory = false;
        bytes_used = 0;
    }

    /**
   * Moves the allocator "stack pointer" back to zero,
   * effectively "reclaiming" all allocated memory.
   */
    void clear() noexcept
    {
#if DEBUG
        std::fill (raw_data.begin(), raw_data.begin() + (std::ptrdiff_t) bytes_used, std::byte { 0xDD });
#endif
        bytes_used = 0;
    }

    /** Returns the number of bytes currently being used */
    [[nodiscard]] size_t get_bytes_used() const noexcept { return bytes_used; }

    /**
   * Allocates a given number of bytes.
   * The returned memory will be un-initialized, so be sure to clear it manually
   * if needed.
   */
    void* allocate_bytes (size_t num_bytes, size_t alignment = 1)
    {
        auto* pointer = snap_pointer_to_alignment (raw_data.data() + bytes_used, alignment);
        const auto bytes_increment = static_cast<size_t> (std::distance (raw_data.data() + bytes_used, pointer + num_bytes));

        if (bytes_used + bytes_increment > raw_data.size())
        {
            assert (false);
            return nullptr;
        }

        bytes_used += bytes_increment;
        return pointer;
    }

    /**
   * Allocates space for some number of objects of type T
   * The returned memory will be un-initialized, so be sure to clear it manually
   * if needed.
   */
    template <typename T, typename IntType>
    T* allocate (IntType num_Ts, size_t alignment = alignof (T))
    {
        return static_cast<T*> (allocate_bytes ((size_t) num_Ts * sizeof (T), alignment));
    }

    /**
   * Returns a span of type T, and size count.
   * The returned memory will be un-initialized, so be sure to clear it manually
   * if needed.
   */
    template <typename T, typename IntType>
    auto make_span (IntType count, size_t alignment = default_byte_alignment)
    {
        return std::span { allocate<T> (count, alignment),
                           static_cast<size_t> (count) };
    }

    /** Returns a pointer to the internal buffer with a given offset in bytes */
    template <typename T, typename IntType>
    T* data (IntType offset_bytes) noexcept
    {
        return reinterpret_cast<T*> (raw_data.data() + offset_bytes);
    }

    /**
   * Creates a "frame" for the allocator.
   * Once the frame goes out of scope, the allocator will be reset
   * to whatever it's state was at the beginning of the frame.
   */
    struct Frame
    {
        Frame() = default;
        explicit Frame (Memory_Arena& allocator)
            : alloc (&allocator), bytes_used_at_start (alloc->bytes_used) {}

        ~Frame()
        {
            if (alloc != nullptr)
                alloc->bytes_used = bytes_used_at_start;
        }

        Memory_Arena* alloc = nullptr;
        size_t bytes_used_at_start = 0;
    };

    /** Creates a frame for this allocator */
    auto create_frame() { return Frame { *this }; }

    void reset_to_frame (const Frame& frame)
    {
        assert (frame.alloc == this);
        bytes_used = frame.bytes_used_at_start;
    }

private:
    void free_owned_memory() noexcept
    {
        if (owns_memory && raw_data.data() != nullptr)
            chowdsp::fft::aligned_free (raw_data.data());
        raw_data = {};
        owns_memory = false;
    }

    std::span<std::byte> raw_data {};
    size_t bytes_used = 0;
    bool owns_memory = false;
};
} // namespace hpss
