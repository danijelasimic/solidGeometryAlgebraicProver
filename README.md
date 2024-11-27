### Constructions

---

- `A = make_point`
- `M = make_midpoint A B`
- `A = make_point_on_line l`
- `A = make_point_on_line P Q`
- `A = make_line_intersection l1 l2`
- `A = make_line_intersection P Q R S`
- `A = make_point_ratio P Q 1 2`  
  *(This means that `AP/AQ = 1/2`. The point `A` is between `P` and `Q`.)*
- `A = make_point_ratio P Q m n`
- `A = make_point_ratio P Q -1 2`  
  *(This means that `AP/AQ = 1/2`. The point `P` is between `A` and `Q`. This is the only command where negative numbers are used.)*
- `A = make_point_in_plane pi`  
  *(Arbitrary point in plane `pi`.)*

---

- `A = translate B 20 30 5`
- `A = translate B a b c`

---

- `A = make_line_plane_intersection l pi`
- `A = make_line_plane_intersection M N pi`
- `A = make_line_plane_intersection l P Q R`
- `A = make_line_plane_intersection M N P Q R`

---

- `A = make_foot_on_plane S pi`  
  *(Projection of point `S` onto plane `pi`.)*
- `Aproj = make_point_projection A pi`  
  *(Projection of point `A` onto plane `pi`.)*  
  Formula: `Aproj = A − ((X−P)⋅n / (n⋅n))⋅n`

---

- `l = make_line A B`
- `l = make_line_orthogonal_on_plane A pi`
- `l = make_line_orthogonal_on_plane A P Q R`

---

- `A B C D A1 B1 C1 D1 = make_cube`
- `A B C D A1 B1 C1 D1 = make_cube 2`  
  *(This means the side length is `2`.)*
- `a = make_cube`

---

- `A B C D = make_tetrahedron`
- `a = make_tetrahedron`
- `A B C D = make_regular_tetrahedron`
- `a = make_regular_tetrahedron`

---

- `A1 A2 A3 A4 A5 A6 = make_regular_hexagon`
- `a = make_regular_hexagon`

---

- `A B C D = make_parallelogram`

---

- `s = make_unit_sphere`  
  *(Creates a unit sphere in canonical position.)*
- `s = make_sphere`
- `s = make_sphere O r`
- `s = make_sphere O number`
- `s = make_sphere A B C D`
- `s O r = make_sphere A B C D`

---

- `n = distance A B`
- `n = square_distance A B`
- `n = 56`
- `n = make_number`
- `n = 56 * n2`  
  *(Number multiplied by an existing number.)*
- `n = m * n`
- `n = a + b`
- `n = n1 / n2`
- `n = distance M A B`  
  *(Distance from point `M` to the line determined by `A` and `B`.)*

---

- `pi = make_plane`  
  *(Creates an arbitrary plane.)*
- `pi = make_plane A B C`
- `pi = make_plane_orthogonal_on_plane tetha l`
- `pi = make_plane_orthogonal_on_plane tetha A B`
- `pi = make_plane_orthogonal_on_line A l`
- `pi = make_plane_tangent_on_sphere s A`

### Assertions

---

- `point_on_line A l`
- `point_on_line A P Q`

---

- `line_intersection A l1 l2`
- `line_intersection A P Q R S`

---

- `orthogonal_lines A B C D`
- `orthogonal_lines l1 l2`

---

- `parallel_lines A B C D`
- `parallel_lines l1 l2`

---

- `equal_points A B`
- `midpoint M A B`
- `point_segment_ratio M C B 20 20`

---

- `congruent A B C D`
- `segments_in_ratio A B C D 2 3`
- `distance A B = distance C D`

---

- `point_in_plane A pi`
- `point_in_plane A P Q R`

---

- `parallel_planes pi teta`
- `parallel_planes A B C P Q R`
- `parallel_planes pi P Q R`

---

- `orthogonal_planes pi teta`
- `orthogonal_planes A B C P Q R`
- `orthogonal_planes pi P Q R`

---

- `parallel_line_plane l pi`
- `parallel_line_plane l P Q R`
- `parallel_line_plane A B pi`
- `parallel_line_plane A B P Q R`

---

- `orthogonal_line_plane l pi`
- `orthogonal_line_plane l P Q R`
- `orthogonal_line_plane A B pi`
- `orthogonal_line_plane A B P Q R`

---

- `line_in_plane l pi`
- `line_in_plane A B pi`
- `line_in_plane l F E G`
- `line_in_plane A B P Q R`

---

- `not_skew p q`
- `not_skew A B P Q`

---

- `point_on_sphere A s`

---

- `line_plane_intersection A l pi`
- `line_plane_intersection A M N pi`
- `line_plane_intersection A l P Q R`
- `line_plane_intersection A M N P Q R`

---

- `equal_angles A O1 B C O2 D`

---

- `equal_numbers n1 n2`

---

- `collinear A B C`
