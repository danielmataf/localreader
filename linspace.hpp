
template<typename T>
auto linspace(T start, T end, const int num) -> std::vector<T> {
    std::vector<T> v;
    if (num == 0) { return v; }
    if (num == 1) {
        v.push_back(start);
        return v;
    }

    double delta = static_cast<double>(end - start) / (num - 1);

    for (int i = 0; i < num - 1; ++i) {
        v.push_back(start + delta * i);
    }
    v.push_back(end);
    return v;
}
